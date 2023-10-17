#ifndef WAVELETLIB_H
#define WAVELETLIB_H

#include "waveletLib_global.h"
#include<QVector>
class WAVELETLIB_EXPORT WaveletLib
{
public:
    WaveletLib();
    double morletWavelet(double t);
    QVector<double> continuousWaveletTransform(const QVector<double>& signal, const QVector<double>& scales, const QVector<double>& positions);
    QVector<double> getRandomNoiseSignal(double amplitude, int sampleNum);
    QVector<double> getFourierSineSignal(double amplitude, double amplitudeFourier, double frequency, int sampleNum, int cycles);
    QVector<double> getFourierCoSineSignal(double amplitude, double amplitudeFourier, double frequency, int sampleNum, int cycles);
    QVector<double> getSquareSignal(double amplitude, double frequency, int sampleNum, int cycles);
    QVector<double> getSineSignal(double amplitude, double frequency, int sampleNum);
    QVector<double> getCoSineSignal(double amplitude, double frequency, int sampleNum);
    QVector<double> getToothSawSignal(double amplitude, double frequency, int sampleNum, int cycles);
    QVector<double> performDFT(QVector<double> signal);
    QVector<QVector<double>> performDWT(QVector<double> signal,QVector<double> phi, QVector<double> psi);
};

#endif // WAVELETLIB_H
