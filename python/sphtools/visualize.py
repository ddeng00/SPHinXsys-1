import argparse
import signal
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
from matplotlib import cm
from matplotlib.colors import ListedColormap
from PyQt5 import Qt
from pyvistaqt import QtInteractor

from .utils import ScreenExporter, _natural_sort_path_key

pv.set_plot_theme('document')
pv.global_theme.multi_samples = 8
PNG_EXPORT_BASE_DIR = Path('./sphviz')


class MainWindow(Qt.QMainWindow):
    point_size = 15
    frame_delay_msec = 20
    default_input_scalar = ' Voltage '
    ScalarRange = {
        'Voltage': (-80, 20),  # (0, 1),
        'GateVariable': (0, 1),
        'ActiveContractionStress': (0, 0.01)
    }
    grays = cm.get_cmap("gray_r")
    cmap = ListedColormap(grays(np.linspace(0, 0.7, 255)), 'MyGrays')
    # cmap = 'bone' # 'gray_r' # 'bone_r'

    # Opacity of the mesh. If a single float value is given, it will be the global opacity of the mesh and uniformly
    # applied everywhere - should be between 0 and 1. A string can also be specified to map the scalars range to a
    # predefined opacity transfer function (options include: 'linear', 'linear_r', 'geom', 'geom_r').
    opacity = None

    GLYPH_ELEM_SKIP = 400
    GLYPH_SCALE_FACTOR = 50

    def __init__(self, path=None, parent=None, show=False, make_pngs=False, make_video=False):
        super().__init__(parent)
        self.title = 'SPH Viewer'
        self._create_gui(parent)

        self.t = 0
        self.timer = Qt.QTimer()
        self.timer.timeout.connect(self.load_next_frame)
        self.timer.setInterval(self.frame_delay_msec)

        self.files = None
        self.tlabel = None

        self.point_cloud = None
        self.actor_pc = None

        self.point_cloud2 = None
        self.actor_pc2 = None
        self.show_pc2 = True

        self.actor_glyth = None
        self.show_glyth = not self.show_pc2

        self.available_input_scalars = []
        self.input_scalar = self.default_input_scalar
        self.scalar = self.input_scalar.strip()

        self.screen_exporter = ScreenExporter()

        if path:
            self.load_dir(path)
        if make_pngs:
            self.toggle_png_sequence()
        if make_video:
            self.toggle_video_recording()

        # Add more lights
        # [print(l) for l in self.plotter.renderer.lights]
        diffuse_color = (0.9998, 0.9998, 0.9998)
        specular_color = (0.9998, 0.9998, 0.9998)
        light = pv.Light(
            light_type='scene light',
            position=(0, -100, -100),
            intensity=0.214286)
        light.diffuse_color = diffuse_color
        light.specular_color = specular_color
        self.plotter.add_light(light)
        light = light.copy()
        light.position = (0, -100, +100)
        self.plotter.add_light(light)

        self.plotter.camera_position = [
            (-87.78549744518088, 68.02895397797286, -356.50205823239264),
            (-10.460600852966309, -34.98706881701946, 0.07710075378417969),
            (0.06570522273760472, 0.9623441385801387, 0.26377373380504077)
        ]

        # Add glyph_widgets
        self.plotter.subplot(0, 1)
        self.plotter.add_checkbox_button_widget(self.toogleGlyph, size=25)

        self.setWindowTitle(self.title)
        w, h = 1386, 742
        self.resize(w, h)
        self.show()
        if not show:  # not calling show() will ignore the window resize() call
            self.hide()

    def _create_gui(self, parent=None):
        Qt.QMainWindow.__init__(self, parent)
        # create the frame
        self.frame = Qt.QFrame()
        vlayout = Qt.QVBoxLayout()

        # add the pyvista interactor object
        # noinspection PyTypeChecker
        self.plotter = QtInteractor(self.frame, shape=(
            1, 2), lighting="light kit")  # "none"
        # Force same view for both subplots
        self.plotter.renderers[1].camera = self.plotter.renderers[0].camera
        # self.plotter.renderer.enable_eye_dome_lighting()
        # camerastyle = vtk.vtkInteractorStyleTrackballCamera()
        # self.plotter.iren.SetInteractorStyle(camerastyle) #https://vtk.org/Wiki/VTK/Examples/Python/InteractorStyleTrackballCamera
        vlayout.addWidget(self.plotter.interactor)

        self.frame.setLayout(vlayout)
        self.setCentralWidget(self.frame)

        # simple menu to demo functions
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        loadButton = Qt.QAction('Load Directory', self)
        loadButton.setShortcut('Ctrl+O')
        loadButton.triggered.connect(self.load_dir_dialog)
        fileMenu.addAction(loadButton)
        movieButton = Qt.QAction('Create Movie', self)
        movieButton.setShortcut('Ctrl+M')
        movieButton.triggered.connect(self.toggle_video_recording)
        fileMenu.addAction(movieButton)
        screenshotButton = Qt.QAction('Toogle PNG Screenshots', self)
        screenshotButton.setShortcut('Ctrl+P')
        screenshotButton.triggered.connect(self.toggle_png_sequence)
        fileMenu.addAction(screenshotButton)
        cameraPosButton = Qt.QAction('Print Camera Position', self)
        cameraPosButton.setShortcut('Ctrl+C')
        cameraPosButton.triggered.connect(self.print_camera_position)
        fileMenu.addAction(cameraPosButton)
        toogleGlyphButton = Qt.QAction('Toogle Glyphs', self)
        toogleGlyphButton.setShortcut('Ctrl+G')
        toogleGlyphButton.triggered.connect(self.toogleGlyph)
        fileMenu.addAction(toogleGlyphButton)
        exitButton = Qt.QAction('Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        # noinspection PyTypeChecker
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)

        toolbar = self.addToolBar("Main")
        self.playpause = Qt.QPushButton("Pause")
        self.playpause.setShortcut("Space")
        toolbar.addWidget(self.playpause)
        self.playpause.clicked.connect(self.toggle_playpause)

        self.comboBox = Qt.QComboBox()
        self.comboBox.activated.connect(self.set_active_scalar)
        toolbar.addWidget(self.comboBox)

        frame_delay_slider = Qt.QSlider(Qt.Qt.Horizontal)
        frame_delay_slider.setValue(self.frame_delay_msec)
        frame_delay_slider.setMinimum(0)
        frame_delay_slider.setMaximum(100)
        frame_delay_slider.valueChanged.connect(self.set_frame_delay)
        toolbar.addWidget(frame_delay_slider)

        pointsize_slider = Qt.QSlider(Qt.Qt.Horizontal)
        pointsize_slider.setValue(self.point_size)
        pointsize_slider.setMinimum(1)
        pointsize_slider.setMaximum(40)
        pointsize_slider.valueChanged.connect(self.set_point_size)
        toolbar.addWidget(pointsize_slider)

        self.time_slider = Qt.QSlider(Qt.Qt.Horizontal)
        self.time_slider.setMinimum(0)
        self.time_slider.setTickInterval(0)
        self.time_slider.valueChanged.connect(self.set_t)
        toolbar.addWidget(self.time_slider)

    def set_glyph_scale(self, val):
        self.GLYPH_SCALE_FACTOR = val

    def set_glyph_skip(self, val):
        self.GLYPH_ELEM_SKIP = int(val)

    def toggle_playpause(self):
        if self.timer.isActive():
            self.timer.stop()
            self.playpause.setText("Play")
        elif self.actor_pc:
            self.timer.setInterval(self.frame_delay_msec)
            self.timer.start()

    def set_t(self):
        self.t = self.time_slider.value()
        self.load(self.t)

    def set_frame_delay(self, val):
        self.frame_delay_msec = val
        self.timer.setInterval(val)

    def set_point_size(self, val):
        self.point_size = val
        if self.actor_pc:
            self.actor_pc.GetProperty().SetPointSize(val)
        if self.actor_pc2:
            self.actor_pc2.GetProperty().SetPointSize(val)

    def set_active_scalar(self, val):
        self.input_scalar = self.available_input_scalars[val]
        self.scalar = self.input_scalar.strip()
        self.point_cloud = pv.PolyData(self.point_cloud.points)
        self.load(self.t)
        self.reload_mesh()

    def load_dir_dialog(self):
        self.timer.stop()
        name = Qt.QFileDialog.getExistingDirectory(self, "Select Directory")
        if name:
            self.t = 0
            self.load_dir(Path(name))
        elif self.actor_pc:
            self.timer.start()

    def locate_files(self, path: Path):
        self.files = sorted(
            path.glob('*.npz'), key=_natural_sort_path_key)
        if not self.files:
            msg = Qt.QMessageBox()
            msg.setIcon(Qt.QMessageBox.Critical)
            msg.setText("Unable to find any suitable files")
            msg.setWindowTitle("Error")
            msg.exec_()

    def reload_mesh(self):
        if self.actor_pc:
            self.plotter.remove_actor(self.actor_pc)
        if self.actor_pc2:
            self.plotter.remove_actor(self.actor_pc2)
            self.actor_pcw2 = None
        if self.actor_glyth:
            self.plotter.remove_actor(self.actor_glyth)
            self.actor_glyth = None

        self.plotter.subplot(0, 0)
        self.actor_pc = self.plotter.add_mesh(self.point_cloud,
                                              render_points_as_spheres=True,
                                              point_size=self.point_size,
                                              lighting=True,
                                              opacity=self.opacity,
                                              cmap=self.cmap,
                                              clim=self.ScalarRange.get(
                                                  self.scalar, (0, 0)),
                                              scalar_bar_args=dict(
                                                  fmt="%.1f",
                                                  color="black"
                                                  # interactive=True
                                              ))

        # self.plotter.scalar_bar.SetAnnotationLeaderPadding(16)
        # self.plotter.scalar_bar.SetTitle("Membrane Potential [mV]")
        # self.tlabel = self.plotter.add_text("test", position=(0.05, 0.05), font_size=11, viewport=True)

        self.plotter.subplot(0, 1)
        cw = 0.95
        if self.show_pc2:
            self.actor_pc2 = self.plotter.add_mesh(self.point_cloud2,
                                                   render_points_as_spheres=True,
                                                   point_size=self.point_size,
                                                   lighting=True,
                                                   opacity=self.opacity,
                                                   # "white"
                                                   color=(cw, cw, cw),
                                                   clim=[0, 0])

    def load_dir(self, path: Path):
        print(f"Loading {path} ...")
        self.locate_files(path)

        d = np.load(self.files[0], mmap_mode='r')
        self.point_cloud = pv.PolyData(d['XYZ'])
        # noinspection PyUnresolvedReferences
        self.available_input_scalars = ['None'] + d.files[1:]

        if self.input_scalar not in self.available_input_scalars:
            if self.input_scalar.strip() in self.available_input_scalars:
                self.input_scalar = self.input_scalar.strip()
            else:
                self.input_scalar = self.available_input_scalars[0]
        self.comboBox.clear()
        self.comboBox.addItems(self.available_input_scalars)
        self.comboBox.setCurrentIndex(
            self.available_input_scalars.index(self.input_scalar))

        self.point_cloud2 = pv.PolyData(self.point_cloud.points)

        self.time_slider.setMaximum(len(self.files))

        if self.t >= len(self.files):
            self.t = 0
        self.load(self.t)
        self.reload_mesh()
        self.screen_exporter.set_params(
            PNG_EXPORT_BASE_DIR / self.files[0].parent.name, len(self.files))

        self.timer.start()

    def load_next_frame(self):
        self.t += 1
        if (self.t < 0) or (self.t >= len(self.files)):
            self.t = 0
        self.load(self.t)

        self.screen_exporter.tick(self.t, self.plotter)

    # noinspection PyTypeChecker
    def load_file(self, path: Path, input_scalar: str):
        d = np.load(path, mmap_mode='r')
        if input_scalar == 'None':
            return d['XYZ'], None
        else:
            return d['XYZ'], d[input_scalar]

    def load(self, t):
        if (t < 0) or (t >= len(self.files)):
            t = 0

        new_points, new_scalar = self.load_file(
            self.files[t], self.input_scalar)
        self.point_cloud.points = new_points
        if self.show_pc2 is not None:
            self.point_cloud2.points = self.point_cloud.points

        if self.input_scalar != 'None':
            if self.scalar == "Voltage":
                self.point_cloud.point_arrays[self.scalar] = 100 * \
                    new_scalar - 80
            else:
                self.point_cloud.point_arrays[self.scalar] = new_scalar

        if self.show_glyth:
            if t + 1 < len(self.files) and t > 0:
                points_t1, _ = self.load_file(self.files[t + 1], 'None')
                vectors = points_t1 - new_points

                pc = pv.PolyData(new_points[::self.GLYPH_ELEM_SKIP])
                pc.vectors = vectors[::self.GLYPH_ELEM_SKIP]
                glyph = pc.glyph(factor=self.GLYPH_SCALE_FACTOR)
                if self.actor_glyth is not None:
                    self.actor_glyth.GetMapper().SetInputData(glyph)
                else:
                    self.plotter.subplot(0, 1)
                    self.actor_glyth = self.plotter.add_mesh(
                        glyph, color='red')
        else:
            self.plotter.remove_actor(self.actor_glyth)
            self.actor_glyth = None

        if self.tlabel is not None:
            self.tlabel.SetInput(f"{t * 12.9:.0f} ms")

        self.time_slider.blockSignals(True)
        self.time_slider.setValue(t)
        self.time_slider.blockSignals(False)

    def toggle_video_recording(self):
        self.screen_exporter.toogle_movie(self.plotter)
        if self.screen_exporter.make_movie:
            self.t = -1

    def toggle_png_sequence(self):
        self.screen_exporter.toogle_png_sequence()
        if self.screen_exporter.make_pngs:
            self.t = -1

    def toogleGlyph(self, _):
        self.show_glyth = not self.show_glyth
        self.show_pc2 = not self.show_pc2
        self.reload_mesh()

        if self.show_glyth:
            self.plotter.subplot(0, 1)
            self.plotter.add_slider_widget(self.set_glyph_scale, (0, 200),
                                           value=self.GLYPH_SCALE_FACTOR, title="Glyph Scale",
                                           pointa=(.1, .08), pointb=(.5, .08),
                                           style='modern')
            self.plotter.add_slider_widget(self.set_glyph_skip, (1, 1000),
                                           value=self.GLYPH_ELEM_SKIP, title="Glyph Skip",
                                           pointa=(.55, .08), pointb=(.95, .08),
                                           style='modern', fmt="%.0f")
            self.plotter.subplot(0, 0)
        else:
            self.plotter.clear_slider_widgets()
            # for widget in self.glyph_widgets: widget.Off()
            # self.glyph_widgets = []

    def print_camera_position(self):
        print(self.plotter.camera_position)


def dir_path_validator(value):
    p = Path(value)
    if not p.exists() or not p.is_dir():
        msg = "{value} is not valid path to an existing directory".format(
            value=value)
        raise argparse.ArgumentTypeError(msg)
    return p.resolve()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'dir', nargs='?', type=dir_path_validator, default=None)
    parser.add_argument('--hide', action='store_true')
    parser.add_argument('--pngs', action='store_true')
    parser.add_argument('--video', action='store_true')
    args, unparsed_args = parser.parse_known_args()

    app = Qt.QApplication(sys.argv[:1] + unparsed_args)
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    # noinspection PyUnusedLocal
    window = MainWindow(path=args.dir, show=not args.hide,
                        make_pngs=args.pngs, make_video=args.video)
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
