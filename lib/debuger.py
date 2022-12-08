import cProfile, io, pstats
import os

PATH = 'debug'

def profile(func):

    if(not os.path.exists(PATH)):
        os.makedirs(PATH)
    
    def wrapper(*args, **kwargs):
        func_name = func.__name__ + ".pfl"
        
        prof = cProfile.Profile()
        retval = prof.runcall(func, *args, **kwargs)
        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(prof, stream=s).strip_dirs().sort_stats(sortby)
        ps.print_stats()
        
        filename = os.path.join(PATH, func_name)
        with open(filename, "w") as f:
            f.write(s.getvalue())
        return retval

    return wrapper