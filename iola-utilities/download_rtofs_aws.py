#!/usr/bin/env python3
import os
import argparse
import concurrent.futures
import subprocess
from pathlib import Path

RTOFS_LIST = """
rtofs_glo.t00z.n-01.archs.a.tgz
rtofs_glo.t00z.n-01.archs.b
rtofs_glo.t00z.n-02.archs.a.tgz
rtofs_glo.t00z.n-02.archs.b
rtofs_glo.t00z.n-03.archs.a.tgz
rtofs_glo.t00z.n-03.archs.b
rtofs_glo.t00z.n-04.archs.a.tgz
rtofs_glo.t00z.n-04.archs.b
rtofs_glo.t00z.n-05.archs.a.tgz
rtofs_glo.t00z.n-05.archs.b
rtofs_glo.t00z.n-06.archs.a.tgz
rtofs_glo.t00z.n-06.archs.b
rtofs_glo.t00z.n-06.archv.a.tgz
rtofs_glo.t00z.n-06.archv.b
rtofs_glo.t00z.n-07.archs.a.tgz
rtofs_glo.t00z.n-07.archs.b
rtofs_glo.t00z.n-08.archs.a.tgz
rtofs_glo.t00z.n-08.archs.b
rtofs_glo.t00z.n-09.archs.a.tgz
rtofs_glo.t00z.n-09.archs.b
rtofs_glo.t00z.n-10.archs.a.tgz
rtofs_glo.t00z.n-10.archs.b
rtofs_glo.t00z.n-11.archs.a.tgz
rtofs_glo.t00z.n-11.archs.b
rtofs_glo.t00z.n-12.archs.a.tgz
rtofs_glo.t00z.n-12.archs.b
rtofs_glo.t00z.n-12.archv.a.tgz
rtofs_glo.t00z.n-12.archv.b
rtofs_glo.t00z.n-13.archs.a.tgz
rtofs_glo.t00z.n-13.archs.b
rtofs_glo.t00z.n-14.archs.a.tgz
rtofs_glo.t00z.n-14.archs.b
rtofs_glo.t00z.n-15.archs.a.tgz
rtofs_glo.t00z.n-15.archs.b
rtofs_glo.t00z.n-16.archs.a.tgz
rtofs_glo.t00z.n-16.archs.b
rtofs_glo.t00z.n-17.archs.a.tgz
rtofs_glo.t00z.n-17.archs.b
rtofs_glo.t00z.n-18.archs.a.tgz
rtofs_glo.t00z.n-18.archs.b
rtofs_glo.t00z.n-18.archv.a.tgz
rtofs_glo.t00z.n-18.archv.b
rtofs_glo.t00z.n-19.archs.a.tgz
rtofs_glo.t00z.n-19.archs.b
rtofs_glo.t00z.n-20.archs.a.tgz
rtofs_glo.t00z.n-20.archs.b
rtofs_glo.t00z.n-21.archs.a.tgz
rtofs_glo.t00z.n-21.archs.b
rtofs_glo.t00z.n-22.archs.a.tgz
rtofs_glo.t00z.n-22.archs.b
rtofs_glo.t00z.n-23.archs.a.tgz
rtofs_glo.t00z.n-23.archs.b
rtofs_glo.t00z.n-24.archs.a.tgz
rtofs_glo.t00z.n-24.archs.b
rtofs_glo.t00z.n-24.archv.a.tgz
rtofs_glo.t00z.n-24.archv.b
rtofs_glo.t00z.n00.archs.a.tgz
rtofs_glo.t00z.n00.archs.b
rtofs_glo.t00z.n00.archv.a.tgz
rtofs_glo.t00z.n00.archv.b
rtofs_glo.t00z.n00.restart.a.tgz
rtofs_glo.t00z.n00.restart.b
""".strip().splitlines()


def download_file(args):
    url, filename = args
    print(f"Downloading: {filename}")
    cmd = ["wget", "-q", url]  # -q for quieter output
    result = subprocess.run(cmd)
    return filename, result.returncode

def extract(tgz):
    print(f"Extracting {tgz}")
    subprocess.run(["tar", "-xzvf", str(tgz)])
    subprocess.run(["rm", str(tgz)])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("date_input", help="YYYYMMDD")
    args = parser.parse_args()

    date_input = args.date_input
    Path(date_input).mkdir(parents=True, exist_ok=True)

    base_url = f"https://noaa-nws-rtofs-pds.s3.amazonaws.com/rtofs.{date_input}"
    print("Base URL:", base_url)

    # Build download tasks
    tasks = [(f"{base_url}/{fname}", fname) for fname in RTOFS_LIST]

    # Download in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=48) as executor:
        results = executor.map(download_file, tasks)

    # Move downloaded files
    for fname, ret in results:
        if ret == 0 and os.path.exists(fname):
            os.rename(fname, f"rtofs.{date_input}/{fname}")
        else:
            print(f"Failed: {fname}")

    # Untar inside the folder
    os.chdir(date_input)
    tgzs=[tgz for tgz in Path(".").glob("*.tgz")]
    # Extract in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=48) as executor:
        results = executor.map(extract, tgzs)
    print("Done.")

if __name__ == "__main__":
    main()

