# NAMD.coor-Amber.parm-to-dms
You can use these files to convert a system simulated in NAMD with Amber parameters into an Anton/Desmond- compatible .dms file.
see note1 to convert the namd simulation output into an amber simulation output, 
then see note2 to take the new amber output file with the original amber parameterizations (.top file), and to convert them into .cms format.
then the mae2dms function in anton will convert .cms to anton-readable .dms
