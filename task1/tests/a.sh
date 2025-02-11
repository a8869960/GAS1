echo "C_rho = mu ="
for ((tau = 0.1; tau <= 0.0001; tau /= 10)) do
for ((h = 0.1; h <= 0.0001; h /= 10))
do ./a 0.1 1 0.1 0.1 l a; done; done
