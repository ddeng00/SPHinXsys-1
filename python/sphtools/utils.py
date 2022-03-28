import re
import argparse
from pathlib import Path


def _natural_sort_path_key(path: Path, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in _nsre.split(path.name)]

def dir_path_validator(value):
    p = Path(value)
    if not p.exists() or not p.is_dir():
        msg = "{value} is not valid path to an existing directory".format(
            value=value)
        raise argparse.ArgumentTypeError(msg)
    return p.resolve()

class ScreenExporter():
    def __init__(self, video_fps=45):
        self.make_pngs = False
        self.make_movie = False

        self.export_path = None
        self.max_t = None
        self.video_fps = video_fps

    def start_movie(self, plotter):
        self.make_movie = True
        p = self.export_path.with_suffix('.mp4')
        print(f"creating {p} ...")
        plotter.open_movie(str(p), framerate=self.video_fps)

        if not self.export_path.exists():
            self.export_path.mkdir(exist_ok=True, parents=True)

    def start_png_sequence(self):
        self.make_pngs = True
        if not self.export_path.exists():
            self.export_path.mkdir(exist_ok=True, parents=True)

    def stop(self, plotter=None):
        self.make_pngs = False
        if self.make_movie:
            if plotter is None:
                raise RuntimeError("plotter is None")
            plotter.mwriter.close()
            self.make_movie = False
            print(f"Video export finished")

    def tick(self, t, plotter):
        if self.make_pngs:
            self._screenshot(self.export_path / f"{t:04d}.png", plotter)
            if t == self.max_t - 1:
                self.stop()

        if self.make_movie:
            plotter.write_frame()
            if t == self.max_t - 1:
                self.stop(plotter)

    def toogle_png_sequence(self):
        if self.make_pngs:
            self.stop()
        else:
            self.start_png_sequence()

    def toogle_movie(self, plotter):
        if self.make_movie:
            self.stop()
        else:
            self.start_movie(plotter)

    def set_params(self, export_dir: Path, size: int):
        self.export_path = export_dir
        self.max_t = size

    def _screenshot(self, png_file: Path, plotter):
        plotter.screenshot(png_file, transparent_background=False)
        print(f"Saved screenshot to {png_file}")
