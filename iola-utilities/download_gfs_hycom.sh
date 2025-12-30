#!/bin/bash

# Usage:
#   ./download_gfs.sh 20251126 20251128
#
# This will download GFS 0.25° pgrb2 files for:
#   2025-11-26 00/06/12/18
#   2025-11-27 00/06/12/18
#   2025-11-28 00/06/12/18

start_date=$1
end_date=$2

if [ -z "$start_date" ] || [ -z "$end_date" ]; then
  echo "Usage: $0 YYYYMMDD YYYYMMDD"
  exit 1
fi

# Convert to seconds since epoch for iteration
current=$(date -d "$start_date" +%s)
end=$(date -d "$end_date" +%s)

# Forecast hours you want
cycles=("00" "06" "12" "18")
fhrs=("000" "003" "006")

while [ "$current" -le "$end" ]; do
    ymd=$(date -u -d "@$current" +%Y%m%d)

    for cyc in "${cycles[@]}"; do

        dirname="gfs.${ymd}${cyc}"
        mkdir -p "$dirname"
        cd "$dirname"

        echo "=== Downloading GFS for ${ymd} ${cyc}Z ==="

        for fh in "${fhrs[@]}"; do
            url="https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.${ymd}/${cyc}/atmos/gfs.t${cyc}z.pgrb2.0p25.f${fh}"
            echo "Downloading: $url"
            wget --no-check-certificate "$url"
        done

        cd ..

    done

    # Move to next day (86400 seconds)
    current=$(( current + 86400 ))
done

