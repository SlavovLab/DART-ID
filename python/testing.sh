# python converter.py ~/git/RTLib/dat/ev.IFN.txt -v -o ~/git/RTLib/dat/ev.IFN_c.txt
# 
# python align.py ~/git/RTLib/dat/ev.IFN.txt --input-converted ~/git/RTLib/dat/ev.IFN_c.txt -v
# python align.py /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC55_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC57_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC61_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC65_hiFDR/evidence.txt -v 
# 
# python update.py /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC55_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC57_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC61_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC65_hiFDR/evidence.txt -v
# 
# python converter.py /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC55_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC57_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC61_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC65_hiFDR/evidence.txt -v -o ~/git/RTLib/dat/ev_nce_c.txt -e ~/git/RTLib/pd_exclude.txt
# 
# python align.py ~/git/RTLib/dat/ev_nce_c.txt -ic
# 
# python align.py ~/git/RTLib/dat/ev_nce_c.txt -ic -v --stan-iters 40000 --rt-distortion 3
# 
# python align.py /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC55_hiFDR/evidence.txt -v --rt-distortion 3 -e ~/git/RTLib/pd_exclude.txt
# 
# python align.py /gd/SingleCell_Data/SCOPE-QE-QC/SQC72/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 3 -f -v --prior-iters 15 --stan-iters 50000 -pf ./alignment_figs
# 
# python align.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 3 -f -v --prior-iters 5 --stan-iters 50000 --remove-exps D[1-4] --filter-pep 0.5 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180420_2
# 
# python converter.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt  -v --remove-exps D[1-4] --filter-pep 0.5 -o /gd/SingleCell_Data/FP17/evidence_c.txt
# 
# python align.py ~/git/RTLib/dat/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 5 -f -v --prior-iters 10 --stan-iters 10000 --filter-pep 0.5 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/SCOPE_20180420_1
# 
# python converter.py ~/git/RTLib/dat/evidence_elite.txt -e ~/git/RTLib/pd_exclude.txt -v --filter-pep 0.5 --filter-num-exps 3 -o ~/git/RTLib/dat/evidence_elite_c.txt
# 
# python align.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 1 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.1 --filter-num-exps 3 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180425_4
# 
# python align.py ~/git/RTLib/dat/evidence_elite_c.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 10 -f -v --prior-iters 15 --stan-iters 100000 --filter-pep 0.5 --filter-num-exps 3 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/SCOPE_20180428_1
# 
# python align.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 1 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.2 --filter-num-exps 3 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180428_6 --type MQ
# 
# python align.py /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC55_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC57_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC61_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC65_hiFDR/evidence.txt -t MQ -f -v --rt-distortion 3 --prior-iters 15 --stan-iters 100000 --filter-pep 0.5 --filter-num-exps 3 --remove-exps SQC55C -o /gd/Slavov_Lab/Albert/RTLib_Alignments/NCE_20180429_2
# 
# python converter.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt -v --filter-pep 0.1 --filter-num-exps 3 -o /gd/SingleCell_Data/FP17/evidence_c.txt
# 
# python align.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 3 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.1 --filter-num-exps 3 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180501_1 --type MQ
# 
# python align.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 0 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.1 --filter-num-exps 3 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180506_mu_lognorm_4 --type MQ
# 
# python align.py /gd/SingleCell_Data/FP17/evidence.txt -e ~/git/RTLib/pd_exclude.txt --rt-distortion 0 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.1 --filter-num-exps 3 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180506_mu_norm_1 --type MQ
# 
# python update.py /gd/SingleCell_Data/FP17/evidence.txt --type MQ -e ~/git/RTLib/pd_exclude.txt --rt-distortion 0 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.1 --filter-num-exps 3 -p /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180506_mu_norm_1 -o /gd/Slavov_Lab/Albert/RTLib_Alignments/FP17_20180506_mu_norm_1 
# 
# python update.py /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC55_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC57_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC61_hiFDR/evidence.txt /gd/Slavov_Lab/SingleCell_Data/SCOPE-QE-QC/SQC65_hiFDR/evidence.txt -t MQ -e ~/git/RTLib/pd_exclude.txt -f -v --rt-distortion 1 --prior-iters 10 --stan-iters 100000 --filter-pep 0.5 --filter-num-exps 3 --remove-exps 55C\|55L\|61A\|61B -o ~/git/RTLib/Alignments/NCE_20180507_3
# 
# python update.py /gd/SingleCell_Data/FP18/evidence.txt --type MQ -e ~/git/RTLib/pd_exclude.txt --rt-distortion 0 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.3 --filter-num-exps 3 -o ~/git/RTLib/Alignments/FP18_20180507_2
# 
# python update.py /gd/SingleCell_Data/FP17/evidence.txt --type MQ -e ~/git/RTLib/pd_exclude.txt --rt-distortion 1 -f -v --prior-iters 10 --stan-iters 100000 --filter-pep 0.5 --filter-num-exps 3 --filter-retention-length 1 -o ~/git/RTLib/Alignments/FP17_20180507_1