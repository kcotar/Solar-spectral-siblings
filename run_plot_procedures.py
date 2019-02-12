from os import system

# py_run = 'solar_sibilings_find_sort_envelopefit.py'
py_run = 'gaia_twins_photometry.py'

arg1 = [5, 7, 10]
# arg2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
arg2 = [None]

for a1 in arg1:
    for a2 in arg2:
        run_str = 'python ' + py_run + ' ' + str(a1)
        if a2 is not None:
            run_str += ' ' + str(a2)
        print 'Running:', run_str
        system(run_str)
