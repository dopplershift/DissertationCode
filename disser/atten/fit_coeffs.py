# Autogenerated on 2013-02-12 05:12:40.254794
from numpy import array
from quantities import cm
assumed_tempc = 10
assumed_wavelength = {'X': array(3.21) * cm, 'S': array(10.0) * cm, 'C': array(5.5) * cm}
assumed_shape = "Brandes"
assumed_canting_width = 10
za_coeffs = {('S', 'V'): array([ 0.00150816,  0.29997535]), ('C', 'H'): array([  3.99174049e-05,   7.70618335e-01]), ('C', 'V'): array([  2.77095504e-05,   8.12085167e-01]), ('S', 'H'): array([ 0.00078871,  0.35900891]), ('X', 'H'): array([ 0.00110675,  0.62143   ]), ('X', 'V'): array([  6.21758502e-04,   6.81336195e-01])}
ka_coeffs = {('S', 'V'): array(0.0193597332371487), ('C', 'H'): array(0.10014759039795704), ('S', 'diff'): array(0.004540782955962982), ('X', 'diff'): array(0.05447737092830401), ('C', 'V'): array(0.07342783867155196), ('S', 'H'): array(0.022597823808004784), ('C', 'diff'): array(0.028855522146278328), ('X', 'H'): array(0.3315601537701545), ('X', 'V'): array(0.27888686163888265)}
sc_coeffs = {('S', 'diff'): array([  8.27820124e+07,  -2.33592404e+00,   6.39257974e+00,
         3.33733415e+00]), ('C', 'H'): array([  2.16095220e-10,   2.01400289e+00,  -2.67685219e+00,
        -1.07009969e+00]), ('X', 'diff'): array([ 0.21760753, -0.28691196,  3.04669013,  1.28819967]), ('C', 'diff'): array([  2.64687508e+03,  -1.31007276e+00,   4.33422913e+00,
         2.31297500e+00]), ('X', 'H'): array([ 0.02494843,  0.27382397, -0.55885854,  0.72062418]), ('S', 'H'): array([  2.08921514e+38,  -8.73679644e+00,   1.57462673e+01,
         9.64223576e+00])}
diff_atten_coeffs = {'X': array(0.15850873597797685), 'S': array(0.17212030958343003), 'C': array(0.2639433784492295)}
del array, cm
