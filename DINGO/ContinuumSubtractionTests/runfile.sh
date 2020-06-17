awk '!/#/' selavy-results_<SB number>_<tile number_<beam number>.txt | awk '{print $8,$9,$19}' >  selavy-results_<SB number>_<tile number_beam<number>_gist.txt
python extractpixspecs.py <SB number> <tile number> <beam number>
python plotfinal.py <SB number> <tile number> <beam number>
pdfjoin finalplots/*<SB number>_<tile number>_<beam number>*png -o finalplots/<SB number>_<tile number>_<beam number>_contsub.pdf
