import cftime
import datetime

__all__ = ['cftime2datetime']
def cftime2datetime(times):
    '''cftime2datetime

    Parameters
    ----------
    times: class in `cftime` or list with cftime arrray.

    Returns
    -------
    out: list with datetime
    '''
    cftime_class = (cftime.Datetime360Day,
                    cftime.DatetimeAllLeap,
                    cftime.DatetimeGregorian,
                    cftime.DatetimeNoLeap,
                    cftime.DatetimeProlepticGregorian,
                    cftime.DatetimeJulian,)

    if isinstance(times, cftime_class):
        times = [times]
    elif isinstance(times, list):
        # check class in times
        if not isinstance(times[0], cftime_class):
            raise Exception("ERROR: list in time is not class of `cftime`.")

    # okay, convert to datetime.
    out = [datetime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in times]
    return out

if __name__ == '__main__':
    import sys, os
    os.chdir(os.path.dirname(__file__))
    t = cftime.DatetimeNoLeap(1980,1,1)
    print(t.year, t.month, t.day, t.hour, t.minute)
