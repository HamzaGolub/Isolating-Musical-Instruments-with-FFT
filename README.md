# Isolating-Musical-Instruments-with-FFT
In this project, musical instruments used in the opening ballad of "I’m Shipping Up To Boston" are isolated using Fast Fourier Transforms and Gabor Transforms.

By using the Gabor Transform, a spectrogram of the frequencies throughout the clip is created. Based on these frequencies, the thresholds of bass versus guitar frequencies can be distinguished. With these thresholds, the Fast Fourier Transform can be run on the clip, after which the frequencies outside of these thresholds can be eliminated. In the case of extracting the baseline/drumbeat, the threshold is between 50 Hz and 200 Hz. The inverse Fourier Transform is used to revert back to the time domain, giving us a soundclip with an isolated baseline/drumbeat. The same process is repeated for the guitar and higher pitched instruments.

MATLAB code written for this project is implemented under [Musical_Freq.m](Musical_Freq.m). The sound clip used for this project can be found in [Clip.mat](Clip.mat). Results are presented in [Isolating_Instruments.pdf](Isolating_Instruments.pdf)
