# mmWave Through Wall Sensing
## Overview
### Motivation
We want to explore the possibility of using mmWave radar sensors to count how many people there are behind the wall.

### Background
The Frequency Modulated Continuous Wave (FMCW) radar is widely utilized for target ranging and velocity estimation. It operates by transmitting a continuous signal that undergoes linear frequency variation over time. This type of signal is called a chirp.

When the transmitted signal hits a target, it reflects back and is received at the sensor end with varying frequencies. The working principle of FMCW radar relies on two factors: the Doppler effect and the delay time of the reflected signal. The Doppler effect causes a frequency shift as a result of relative motion between the radar sensor and the target. By analyzing this frequency shift, the radar system can accurately determine the velocity of the target.

By separating different distance and velocity of target, we can detect the number of target objects behind the screen.

### Implementation
#### Hardware
- xWR6843 mmWave sensor
- DCA1000EVM capture card
- mmWave studio tool (UI for sensor configurations and easy capture)

#### Data capture
- Use the mmWave Studio to capture raw binary data
- Run [rawDataReader.m](rawDataReader.m) to construct raw data in frames, chirps and channels, then performs range FFT to generate radar cube

#### Data Processing
- We explored with different filtering methods to filter out the static white noise
- We experimented with distance / velocity filters, as well as various filtering parameters
- We applied convolution masks to blur data points and decrease over-sampling

#### Algorithm
- [countPeople.m](countPeople.m) uses K-means clustering with Elbow Method to find the number of clusters
- [Euclidean-find-peaks.m](Euclidean-find-peaks.m) uses Euclidean distances to separate the boudaries of maximal clusters

### Results
- Accuracy to count 1-3 people behind wall is 80%
- The application is limited by wall material. Reflective materials such as metal cannot be penetrated with mmWave radar.

## Documentations
- Final result presentation [here](mmWave-through-wall-imaging-presentation.pdf)
- Project report [here](Final_Project_Report.pdf)

## Credits
- Yi-Ting Wu, Akshata Patil, Jesica Gonzalez, Lilian Chan
- Stanford EE292Q Spring 2023
