folder='../../bayesian-ergodic-search/assets/double_integrator3D_TSPN/medium_2D_Behdani_CETSPs/'
args='Results/medium_2D_Behdani_CETSPs/radius_0.25/ 2D 1 14400 V1 BeFS 1 1'
ipython exec/run_instances.py 7 run/exeCVXHULL ${folder} ${args}
ipython exec/run_instances.py 7 run/clarabel_redundant $folder $args
ipython exec/run_instances.py 7 run/clarabel_dropin $folder $args
ipython exec/run_instances.py 7 run/clarabel_recycling $folder $args
ipython exec/run_instances.py 7 run/clarabel_warmstart $folder $args
ipython exec/run_instances.py 7 run/scs_dropin $folder $args
ipython exec/run_instances.py 7 run/scs_recycling $folder $args
ipython exec/run_instances.py 7 run/scs_warmstart $folder $args
ipython exec/run_instances.py 7 run/scs_tridiag_dropin $folder $args
ipython exec/run_instances.py 7 run/scs_tridiag_recycling $folder $args
ipython exec/run_instances.py 7 run/scs_tridiag_warmstart $folder $args
