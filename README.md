# CO685Project
Final project for CO685: The Mathematics of Public-key Cryptography

This is an Implementation in Julia of the Castryck Decru attack on SIKE (https://eprint.iacr.org/2022/975.pdf)

The file curves.jl contains my own implementation of elliptic and hyperelliptic curve arithmetic as well as an implementation of the Weil pairing and Velu's formulas for isogenies which I could not find in Julia

The other files are a direct translation of the sage files of the same name at https://github.com/jack4818/Castryck-Decru-SageMath. 

Line 18 in the file SIKEp434 can be changed to use different parameters.
 
The package Nemo (https://nemocas.github.io/Nemo.jl/latest/) is used for finite field arithmetic 
