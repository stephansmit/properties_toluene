# Properties of Toluene
Scripts to analyse the properties of Toluene

## Requirements
~~~~
python3 -m pip install numpy pandas CoolProp matplotlib
~~~~

## Usage 
Create the module
~~~~
cd fluidmodels_su2
swig -c++ -python3 su2_models.i
python3 setup.py build_ext --inplace
~~~~

Run the script
~~~~
python3 run.py
~~~~

~
