#include "waveletlib.h"
#include<QVector>
#include<math.h>
#include <QRandomGenerator>
#include <complex>
typedef std::complex<double> Complex;

WaveletLib::WaveletLib()
{
}

double WaveletLib::morletWavelet(double t)
{
    double sigma = 0.9;
    double f = 0.5;
    double A = 1.0 / sqrt(sigma * sqrt(2 * M_PI));
    double expPart = exp(-t * t / (2 * sigma * sigma));
    double cosPart = cos(2 * M_PI * f * t);
    double sinPart = sin(2 * M_PI * f * t);
    std::complex<double> morletValue(A * expPart * cosPart, A * expPart * sinPart);
    return abs(morletValue);
//    return std::exp(-t * t / (2 * sigma * sigma)) * std::cos(2 * M_PI * freq * t);
}


QVector<double> WaveletLib::continuousWaveletTransform(const QVector<double>& signal, const QVector<double>& scales, const QVector<double>& positions)
{
    QVector<double> cwt;

    for (double s : scales)
    {
        QVector<double> cwtValues;
        for (double p : positions)
        {
            double sum = 0.0;
            for (int i = 0; i < signal.size(); ++i)
            {
                double t = (i - p) / s;
                sum += signal[i] * morletWavelet(t);
            }
            cwtValues.append(sum);
        }
        cwt.append(cwtValues);
    }

    return cwt;
}


QVector<double> WaveletLib::getRandomNoiseSignal(double amplitude, int sampleNum){
    QVector<double> signal;
    QRandomGenerator rd;
    for(int i =0; i<sampleNum;i++){

        double val;
        if(rd.bounded(0,2) == 0){
            val = rd.generateDouble() * (-amplitude);
        }
        else{
            val = rd.generateDouble() * (amplitude);
        }
        signal.append(val);
    }
    return signal;
}

QVector<double> WaveletLib::performDFT(QVector<double> signal){
    const int numSamples = signal.size();
    int N = signal.size();
    QVector<Complex> spectrum(N);

    for (int k = 0; k < N; ++k) {
        Complex sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            Complex expo = std::polar(1.0, -2 * M_PI * k * n / N);
            sum += signal.at(n) * expo;
        }
        spectrum[k] = sum;
    }

    QVector<double> magnitudeSpectrum;
    for (int i = 0; i < numSamples; ++i) {
        double mag = std::abs(spectrum.at(i));
        magnitudeSpectrum.append(mag);
    }
    return magnitudeSpectrum;
}


QVector<QVector<double>> WaveletLib::performDWT(QVector<double> signal,QVector<double> phi, QVector<double> psi){
    int signalLength = signal.size();
    int filterLength = phi.size();

    QVector<double> dwtPHICoefficients(signalLength);
    QVector<double> dwtPSICoefficients(signalLength);
    QVector<QVector<double>> dwtCoefficients;
    // Convolve the signal with the wavelet filters
    for (int i = 0; i < signalLength; ++i) {
        double approx = 0.0;
        double detail = 0.0;

        for (int j = 0; j < filterLength; ++j) {
            int index = (i - j + filterLength) % filterLength;
            approx += signal[i] * phi[index];
            detail += signal[i] * psi[index];
        }

        dwtPSICoefficients[i] = detail;
        dwtPHICoefficients[i] = approx;
    }
    dwtCoefficients.append(dwtPHICoefficients);
    dwtCoefficients.append(dwtPSICoefficients);
    return dwtCoefficients;
}



QVector<double> WaveletLib::getSquareSignal(double amplitude, double frequency, int sampleNum, int cycles){
    QVector<double> signal;
    const double period = 1.0/frequency;
    for(int c=0; c<cycles;++c){
        for(int i =0; i<sampleNum;i++){
            double time = c * period + i * period / sampleNum;
            //        double val = (t < period / 2.0) ? a : -a;
            double phase = std::fmod(time,period);
            double val = (phase < period / 2.0) ? amplitude : -amplitude;
            signal.append(val);
        }
    }
    return signal;
}

QVector<double> WaveletLib::getFourierSineSignal(double amplitude, double amplitudeFourier, double frequency, int sampleNum, int cycles){
    const double period = 1.0/frequency;
    double offset = amplitude;
    QVector<double> signal;
    for(int c=0; c<cycles;++c){
        for(int i =0; i<sampleNum;i++){
            double val;
            if(i==0 || i == sampleNum/2 || i == sampleNum){
                double t = c * period + i * period / sampleNum;
                double phase = std::fmod(t,period);
                val = (phase < period / 2.0) ? amplitude : -amplitude;
                offset = val;
            }
            else{
                double t = i * 0.01;
                val = amplitudeFourier * std::sin(2.0 * M_PI * frequency * t) + offset;
            }
            signal.append(val);
        }
    }
    return signal;
}

QVector<double> WaveletLib::getFourierCoSineSignal(double amplitude, double amplitudeFourier, double frequency, int sampleNum, int cycles){
    const double period = 1.0/frequency;
    double offset = amplitude;
    QVector<double> signal;
    for(int c=0; c<cycles;++c){
        for(int i =0; i<sampleNum;i++){
            double val;
            if(i==0 || i == sampleNum/2 || i == sampleNum){
                double t = c * period + i * period / sampleNum;
                double phase = std::fmod(t,period);
                val = (phase < period / 2.0) ? amplitude : -amplitude;
                offset = val;
            }
            else{
                double t = i * 0.01;
                val = amplitudeFourier * std::cos(2.0 * M_PI * frequency * t) + offset;
            }
            signal.append(val);
        }
    }
    return signal;
}

QVector<double> WaveletLib::getSineSignal(double amplitude, double frequency, int sampleNum){
    QVector<double> signal;
    for(int i =0; i<sampleNum;i++){
        double t = i * 0.01;
        double val = amplitude * std::sin(2.0 * M_PI * frequency * t);
        signal.append(val);
    }
    return signal;
}

QVector<double> WaveletLib::getCoSineSignal(double amplitude, double frequency, int sampleNum){
    QVector<double> signal;
    for(int i =0; i<sampleNum;i++){
        double t = i * 0.01;
        double val = amplitude * std::cos(2.0 * M_PI * frequency * t);
        signal.append(val);
    }
    return signal;
}

QVector<double> WaveletLib::getToothSawSignal(double amplitude, double frequency, int sampleNum, int cycles){
    QVector<double> signal;
    const double period = 1.0 / frequency;
    const double cycleDuration = period * cycles;

    for(int i =0; i<sampleNum*cycles;++i){
        double t = i *cycleDuration/ (sampleNum * cycles);
        //        double val = (t < period / 2.0) ? a : -a;
        double phase = std::fmod(t,period);
        double val = (phase / period) * 2.0 * amplitude - amplitude;
        signal.append(val);
    }
    return signal;
}
