folder='../../bayesian-ergodic-search/assets/double_integrator3D_TSPN/medium_2D_Behdani_CETSPs/'
args='Results/medium_2D_Behdani_CETSPs/radius_0.25/ 2D 1 14400 V1 BeFS 1 1'
ipython exec/run_instances.py 7 run/exeCVXHULL ${folder} ${args}
ipython exec/run_instances.py 7 run/clarabel_redundant $folder $args
ipython exec/run_instances.py 7 run/clarabel_reduce $folder $args
ipython exec/run_instances.py 7 run/clarabel_reuse $folder $args
ipython exec/run_instances.py 7 run/clarabel_recycle $folder $args
ipython exec/run_instances.py 7 run/scs_reduce $folder $args
ipython exec/run_instances.py 7 run/scs_reuse $folder $args
ipython exec/run_instances.py 7 run/scs_recycle $folder $args
ipython exec/run_instances.py 7 run/scs_tridiag_reduce $folder $args
ipython exec/run_instances.py 7 run/scs_tridiag_reuse $folder $args
ipython exec/run_instances.py 7 run/scs_tridiag_recycle $folder $args
