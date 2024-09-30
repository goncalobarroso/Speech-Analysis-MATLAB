# 🎤 Audio Digit Recognition using MATLAB

This project focuses on the recognition of spoken digits (0-9) in English using signal processing techniques in **MATLAB**. We analyze audio recordings to extract unique temporal and spectral features for effective digit classification.

[![MATLAB](https://img.shields.io/badge/MATLAB-Signal%20Processing-yellow)](https://www.mathworks.com/products/matlab.html)  
[![Dataset](https://img.shields.io/badge/Dataset-AudioMNIST-brightgreen)](https://github.com/soerenab/AudioMNIST)  
[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue)](https://github.com/goncalobarroso/audio-digit-recognition)

## 🚀 Project Overview

The goal is to develop a robust model for identifying spoken digits from audio signals. The **AudioMNIST** dataset was used, consisting of 60 participants, each pronouncing digits from 0 to 9 fifty times.

## 📋 Description

### Features:
- **Temporal Analysis**: Extracting features like energy, amplitude, and signal duration to distinguish digits.
- **Spectral Analysis**: Using **Fast Fourier Transform (FFT)** to analyze the frequency domain and extract the median amplitude spectrum.
- **Short-Time Fourier Transform (STFT)**: Time-varying frequency analysis for capturing complex audio characteristics.

### Dataset:
- **60 participants**, each saying digits **0-9**.
- **48 KHz** sampling rate.
- Each digit is spoken **50 times**, leading to a rich variety of audio data for analysis.
  
## 🛠️ Technologies Used

- **MATLAB** ![MATLAB](https://img.shields.io/badge/-MATLAB-0076A8?logo=mathworks&logoColor=white)  
- **Signal Processing**: Extensive use of MATLAB’s built-in functions for signal analysis.

## 🖼️ Features

- **Temporal Feature Extraction**: Distinguishing digits using their energy, amplitude, and duration.

- **Spectral Analysis**: Using **FFT** to examine dominant frequency components of each digit.

- **STFT Analysis**: Capturing time-frequency data for more detailed recognition.

## 🔬 Methodology

1. **Data Import**: Audio files are imported and organized for efficient processing.
2. **Temporal Analysis**: Investigating features like signal energy and amplitude.
3. **Spectral Analysis**: Using **FFT** for frequency-domain insights.
4. **STFT**: Applying **Short-Time Fourier Transform** for time-varying frequency analysis.
5. **Decision Rules**: Developing logical rules to classify digits based on extracted features.

## 🔗 Useful Links

- **Dataset**: [AudioMNIST on GitHub](https://github.com/soerenab/AudioMNIST)
- **Paper**: [arXiv: Interpreting and Explaining Deep Neural Networks for Classification of Audio Signals](https://arxiv.org/abs/1807.03418)
- **Additional Resources**: [AudioMNIST on Kaggle](https://www.kaggle.com/datasets/sripaadsrinivasan/audio-mnist)

## 📄 License

This project does not use any specific license. Feel free to explore and contribute.
