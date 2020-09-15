awk '!/#/' selavy-image.i.<SB number>.cont.taylor.0.restored.islands.txt | awk '{print $6,$7,$9,$10,$11,$12}' >  selavy-results_<SB number>_gist.txt
python extractpixspecs_top10.py <SB number>
python plotfinal.py <SB number>
