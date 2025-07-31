# cardiac_diffusion
Files needed to generate cardiac diffusion MRI sequences utilizing the MATLAB version of Pulseq as well as example reconstruction code.

The provided MATLAB code requires the open-source pulse sequence framework Pulseq https://pulseq.github.io/ to run.
The file "directions.mat" contains example diffusion directions but can be replaced with any desired q-space sampling scheme.
No guarantees regarding compatability of generated files with any hardware or the suitability for use on volunteers/phantoms/patients are made.
Special Thanks go to Prof. Dr. Maxim Zaitsev and the rest of the Pulseq team for assisting with development. The grappa reconstruction MATLAB file created
by Marh Chiew https://github.com/mchiew/grappa-tools is required to run the reconstruction.

Questions regarding the project can be addressed to eric.arbes@uniklinik-freiburg.de

Please note that the provided code examples are somewhat outdated and feature some hardcoded parameters. It is recommended to use the labeling function in Pulseq to generate data that is easier to handle.
This will be remedied in a future update to the code.

Step by step guide on how to get some initial images:

1. Make sure SequenceGenerationV1, Data Reconstruction V1, directions.m and grappa.m are present in your directory.
2. Leave the parameters in SequenceGenerationV1 unchanged if you don't want to adapt the data reconstruction code and generate three different sequence files by changing the "part" variable in line 65.
3. Run the sequence at your scanner and export the raw data. Note that a ecg trigger signal is required for the sequence to advance, meaning that a synthetic signal should be set at the scanner if phantoms or volunteers without ecg applied are measured.
4. Provide the raw data directory and one of the generated sequences at the start of DataReconstructionV1. Change "numberofaverages" to how many times you ran all sequence parts for averaging.
5. Make sure your raw data directory features the nested subfolders "temp" and "temp2" (This is a legacy step that will be unneeded in the next update)
6. Running the code will yield he b0 and directional images for all slices
