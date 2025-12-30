import os
import pandas as pd

def st(x):
    if x>9:
        return str(x)
    else:
        return '0'+str(x)

start_date='2021-04-08 00:00:00'

start_date=pd.to_datetime(start_date)
end_date=start_date+pd.Timedelta(days=7)

for date in pd.date_range(start_date,end_date,freq='6H'):
    os.system(f'wget --no-check-certificate https://data.rda.ucar.edu/d083003/{date.year}/{date.year}{st(date.month)}/gdas1.fnl0p25.{date.year}{st(date.month)}{st(date.day)}{st(date.hour)}.f00.grib2')
    os.system(f'wget --no-check-certificate https://data.rda.ucar.edu/d083003/{date.year}/{date.year}{st(date.month)}/gdas1.fnl0p25.{date.year}{st(date.month)}{st(date.day)}{st(date.hour)}.f03.grib2')
    
