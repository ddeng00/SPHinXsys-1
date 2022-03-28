import argparse
from pathlib import Path
import numpy as np
from scipy.interpolate import griddata
from multiprocessing import Pool
from .utils import _natural_sort_path_key, dir_path_validator
import time


def dir_path_generator(value):
    p = Path(value)
    if not p.exists():
        p.mkdir(exist_ok=True, parents=True)
    return p.resolve()


def resample(in_dir, out_dir, grid, std_grid, time_step):
    start_time = time.time()
    print(f'Processing {in_dir}')
    files = sorted(in_dir.glob('*.npz'), key=_natural_sort_path_key)
    for i in range(0, len(files), time_step+1):
        # Load data in temporal sequence
        data = [np.load(files[i+j], mmap_mode='r') for j in range(time_step+1)]

        # Interpolate voltage
        xyz0 = data[0]['XYZ']
        voltage = griddata(xyz0, data[0]['Voltage'], grid, fill_value=0)
        # Store voltage in usable matrix
        volt_mat = np.zeros((96, 80, 72))
        volt_mat[std_grid] = voltage.ravel()

        disp_mats = []
        for d in data[1:]:
            # Interpolate displacement
            xyz1 = d['XYZ']
            displacement = griddata(xyz0, np.subtract(
                xyz1, xyz0), grid, fill_value=0)
            # Store displacement in usable matrix
            for j in range(3):
                disp_mat = np.zeros((96, 80, 72))
                disp_mat[std_grid] = displacement[:,j]
                disp_mats.append(disp_mat)
            xyz0 = xyz1

        # Store data
        filename = in_dir.stem + '_' + files[i].stem
        np.savez_compressed(out_dir/filename,
                            Displacement=np.asarray(disp_mats), Voltage=np.asarray([volt_mat]))
    print(f'Finished processing {in_dir} in {time.time()-start_time} seconds')


def main():
    start_time = time.time()

    # Process command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('in_dir', nargs='?',
                        type=dir_path_validator, default=None)
    parser.add_argument('out_dir', nargs='?',
                        type=dir_path_generator, default=None)
    parser.add_argument('time_step', nargs='?', type=int, default=1)
    args, _ = parser.parse_known_args()
    in_dir, out_dir = args.in_dir, args.out_dir
    time_step = args.time_step

    # Generate target coordinates
    xx,yy,zz=np.arange(-55, 35), np.arange(-75, 5), np.arange(-35, 35)
    xx,yy,zz = np.meshgrid(xx,yy,zz)
    grid = xx.ravel(), yy.ravel(), zz.ravel()
    std_grid = np.add(grid[0], 55), np.add(grid[1], 75), np.add(grid[2], 35)

    # Resample generated data
    pool = Pool()
    for path in in_dir.iterdir():
        if path.is_dir():
            pool.apply_async(resample, args=(
                path, out_dir, grid, std_grid, time_step,))
    pool.close()
    pool.join()

    print(f'Total execution time: {time.time() - start_time} seconds')


if __name__ == '__main__':
    main()
