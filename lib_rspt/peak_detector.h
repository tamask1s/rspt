/***************************************************************************
* Copyright 2024 Tamas Levente Kis - tamkis@gmail.com                      *
*                                                                          *
* Licensed under the Apache License, Version 2.0 (the "License");          *
* you may not use this file except in compliance with the License.         *
* You may obtain a copy of the License at                                  *
*                                                                          *
*     http://www.apache.org/licenses/LICENSE-2.0                           *
*                                                                          *
* Unless required by applicable law or agreed to in writing, software      *
* distributed under the License is distributed on an "AS IS" BASIS,        *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
* See the License for the specific language governing permissions and      *
* limitations under the License.                                           *
***************************************************************************/

/**
 * @brief Class for detecting peaks in a signal.
 *
 * This class implements a peak detection algorithm using a combination of
 * bandpass filtering, integration, and adaptive thresholding.
 * It is designed to work with real-time sampled signals, identifying peaks
 * based on a dynamic threshold derived from the signal itself.
 *
 * The detection logic includes:
 * - Bandpass filtering to isolate relevant frequency components.
 * - Integration to smooth out short-term fluctuations.
 * - Adaptive thresholding to dynamically adjust sensitivity.
 *
 * Once a peak is detected, a marker value is returned, otherwise, zero is returned.
 * Peaks are detected with a delay of 220 ms.
 */
class peak_detector
{
    iir_filter_4th_order bandpass_filter_;
    iir_filter_2nd_order integrative_filter_;
    iir_filter_2nd_order threshold_filter_;

    double previous_peak_amplitude_ = 0;
    double previous_sig_val_ = 0;
    bool searching_for_peaks_ = false;
    int samples_after_peak_count_ = 0;
    int sample_indx_ = 0;

    const double sampling_rate_;
    const double marker_val_;
    const double previous_peak_reference_ratio_ = 0.5;
    const double previous_peak_reference_attenuation_ = 25;
    const double peak_attenuation_;
    const double threshold_ratio_ = 1.5;
    const int nr_slope_samples_;

public:
    /**
     * @brief Constructor for the peak detector.
     *
     * Initializes the detector with a given sampling rate and optional marker value.
     * It also configures the internal filters using predefined Butterworth filter settings.
     *
     * @param sampling_rate The sampling rate of the input signal in Hz.
     * @param marker_val The value returned when a peak is detected (default: 1).
     */
    peak_detector(double sampling_rate, double marker_val = 1)
        : sampling_rate_(sampling_rate),
          marker_val_(marker_val),
          peak_attenuation_(1.0 / (1.0 + previous_peak_reference_attenuation_ / sampling_rate)),
          nr_slope_samples_((100.0 * sampling_rate) / 1000.0)
    {
        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 2, sampling_rate, 10, 20);
        create_filter_iir(integrative_filter_.d, integrative_filter_.n, butterworth, low_pass, 2, sampling_rate, 3, 0);
        create_filter_iir(threshold_filter_.d, threshold_filter_.n, butterworth, low_pass, 2, sampling_rate, 0.15, 0);
    }

    /**
     * @brief Processes a new sample and detects peaks.
     *
     * This function applies the filtering and peak detection logic to an incoming sample.
     * The algorithm tracks changes in the signal and determines if a peak has occurred
     * based on an adaptive threshold.
     *
     * @param new_sample The new sample from the signal.
     * @return marker_val_ if a peak is detected, otherwise 0.
     */
    inline double detect(double new_sample, double* peak_sample = 0, double* threshold_sample = 0)
    {
        if (!sample_indx_++)
            bandpass_filter_.init_history_values(new_sample, sampling_rate_);

        double sig_val = bandpass_filter_.filter(new_sample);
        sig_val = integrative_filter_.filter(sig_val * sig_val);
        double threshold = threshold_filter_.filter(sig_val);
        (peak_sample) ? ((*peak_sample) = sig_val) : sample_indx_ = sample_indx_;
        (threshold_sample) ? ((*threshold_sample) = threshold) : sample_indx_ = sample_indx_;

        if (searching_for_peaks_ && (sig_val > threshold * threshold_ratio_) && (previous_sig_val_ > sig_val))
        {
            if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
            {
                previous_peak_amplitude_ = previous_sig_val_;
                samples_after_peak_count_ = 1;
                searching_for_peaks_ = false;
            }
            else
                previous_peak_amplitude_ *= peak_attenuation_;
        }
        else if (previous_sig_val_ < sig_val)
        {
            searching_for_peaks_ = true;
            samples_after_peak_count_ = 0;
        }

        previous_sig_val_ = sig_val;

        if (samples_after_peak_count_)
            ++samples_after_peak_count_;

        if (samples_after_peak_count_ == nr_slope_samples_)
        {
            samples_after_peak_count_ = 0;
            return ((marker_val_ == -1.0) ? sig_val : marker_val_);
        }
        return 0;
    }
};

class peak_detector_1st_order
{
    iir_filter_2nd_order bandpass_filter_;
    iir_filter_1st_order integrative_filter_;
    iir_filter_2nd_order threshold_filter_;

    double previous_peak_amplitude_ = 0;
    double previous_sig_val_ = 0;
    bool searching_for_peaks_ = false;
    int samples_after_peak_count_ = 0;
    int sample_indx_ = 0;

    const double sampling_rate_;
    const double marker_val_;
    const double previous_peak_reference_ratio_ = 0.5;
    const double previous_peak_reference_attenuation_ = 25;
    const double peak_attenuation_;
    const double threshold_ratio_ = 1.5;
    const int nr_slope_samples_;

public:
    /**
     * @brief Constructor for the peak detector.
     *
     * Initializes the detector with a given sampling rate and optional marker value.
     * It also configures the internal filters using predefined Butterworth filter settings.
     *
     * @param sampling_rate The sampling rate of the input signal in Hz.
     * @param marker_val The value returned when a peak is detected (default: 1).
     */
    peak_detector_1st_order(double sampling_rate, double marker_val = 1)
        : sampling_rate_(sampling_rate),
          marker_val_(marker_val),
          peak_attenuation_(1.0 / (1.0 + previous_peak_reference_attenuation_ / sampling_rate)),
          nr_slope_samples_((100.0 * sampling_rate) / 1000.0)
    {
        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 1, sampling_rate, 10, 20);
        create_filter_iir(integrative_filter_.d, integrative_filter_.n, butterworth, low_pass, 1, sampling_rate, 3, 0);
        create_filter_iir(threshold_filter_.d, threshold_filter_.n, butterworth, low_pass, 2, sampling_rate, 0.15, 0);
    }

    /**
     * @brief Processes a new sample and detects peaks.
     *
     * This function applies the filtering and peak detection logic to an incoming sample.
     * The algorithm tracks changes in the signal and determines if a peak has occurred
     * based on an adaptive threshold.
     *
     * @param new_sample The new sample from the signal.
     * @return marker_val_ if a peak is detected, otherwise 0.
     */
    inline double detect(double new_sample, double* peak_sample = 0, double* threshold_sample = 0)
    {
        if (!sample_indx_++)
            bandpass_filter_.init_history_values(new_sample, sampling_rate_);

        double sig_val = bandpass_filter_.filter(new_sample);
        sig_val = integrative_filter_.filter(sig_val * sig_val);
        double threshold = threshold_filter_.filter(sig_val);
        (peak_sample) ? ((*peak_sample) = sig_val) : sample_indx_ = sample_indx_;
        (threshold_sample) ? ((*threshold_sample) = threshold) : sample_indx_ = sample_indx_;

        if (searching_for_peaks_ && (sig_val > threshold * threshold_ratio_) && (previous_sig_val_ > sig_val))
        {
            if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
            {
                previous_peak_amplitude_ = previous_sig_val_;
                samples_after_peak_count_ = 1;
                searching_for_peaks_ = false;
            }
            else
                previous_peak_amplitude_ *= peak_attenuation_;
        }
        else if (previous_sig_val_ < sig_val)
        {
            searching_for_peaks_ = true;
            samples_after_peak_count_ = 0;
        }

        previous_sig_val_ = sig_val;

        if (samples_after_peak_count_)
            ++samples_after_peak_count_;

        if (samples_after_peak_count_ == nr_slope_samples_)
        {
            samples_after_peak_count_ = 0;
            return ((marker_val_ == -1.0) ? sig_val : marker_val_);
        }
        return 0;
    }
};

class peak_detector_offline
{
    iir_filter_2nd_order bandpass_filter_;
    iir_filter_1st_order integrative_filter_;
    iir_filter_1st_order baseline_filter_;
    iir_filter_2nd_order threshold_filter_;

    double previous_peak_amplitude_ = 0;
    double previous_sig_val_ = 0;
    bool searching_for_peaks_ = false;
    int samples_after_peak_count_ = 0;
    int sample_indx_ = 0;

    const double sampling_rate_;
    const double marker_val_;
    const double previous_peak_reference_ratio_ = 0.5;
    const double previous_peak_reference_attenuation_ = 70;
    const double peak_attenuation_;
    const double threshold_ratio_ = 1.5;
    const int nr_slope_samples_;

public:
    /**
     * @brief Constructor for the peak detector.
     *
     * Initializes the detector with a given sampling rate and optional marker value.
     * It also configures the internal filters using predefined Butterworth filter settings.
     *
     * @param sampling_rate The sampling rate of the input signal in Hz.
     * @param marker_val The value returned when a peak is detected (default: 1).
     */
    peak_detector_offline(double sampling_rate, double marker_val = 1)
        : sampling_rate_(sampling_rate),
          marker_val_(marker_val),
          peak_attenuation_(1.0 / (1.0 + previous_peak_reference_attenuation_ / sampling_rate)),
          nr_slope_samples_((100.0 * sampling_rate) / 1000.0)
    {
        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 1, sampling_rate, 15, 25);
        create_filter_iir(integrative_filter_.d, integrative_filter_.n, butterworth, low_pass, 1, sampling_rate, 3, 0);
        create_filter_iir(baseline_filter_.d, baseline_filter_.n, butterworth, low_pass, 1, sampling_rate, 0.5, 0);
        create_filter_iir(threshold_filter_.d, threshold_filter_.n, butterworth, low_pass, 2, sampling_rate, 0.15, 0);
    }

    inline void detect_fw(double* ecg_signal, unsigned int len, double* peak_signal, double* filt_signal, double* threshold_signal)
    {
        bandpass_filter_.init_history_values(ecg_signal[0], sampling_rate_);

        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = bandpass_filter_.filter(ecg_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i] * filt_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            threshold_signal[i] = threshold_filter_.filter(filt_signal[i]);

        for (unsigned int i = 0; i < len; ++i)
        {
            if (searching_for_peaks_ && (filt_signal[i] > threshold_signal[i] * threshold_ratio_) && (previous_sig_val_ > filt_signal[i]))
            {
                if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
                {
                    previous_peak_amplitude_ = previous_sig_val_;
                    samples_after_peak_count_ = 1;
                    searching_for_peaks_ = false;
                }
                else
                    previous_peak_amplitude_ *= peak_attenuation_;
            }
            else if (previous_sig_val_ < filt_signal[i])
            {
                searching_for_peaks_ = true;
                samples_after_peak_count_ = 0;
            }

            previous_sig_val_ = filt_signal[i];

            if (samples_after_peak_count_)
                ++samples_after_peak_count_;

            if (samples_after_peak_count_ == nr_slope_samples_)
            {
                samples_after_peak_count_ = 0;
                peak_signal[i] = ((marker_val_ == -1.0) ? filt_signal[i] : marker_val_);
            }
            else
                peak_signal[i] = 0;
        }
    }

    inline void detect(double* ecg_signal, unsigned int len, double* peak_signal, double* filt_signal, double* threshold_signal, std::vector<unsigned int>* peak_indexes = 0)
    {
        bandpass_filter_.init_history_values(ecg_signal[0], sampling_rate_);
        baseline_filter_.init_history_values(ecg_signal[0], sampling_rate_);

        double* baseline = new double[len];
        for (unsigned int i = 0; i < len; ++i)
            baseline[i] = baseline_filter_.filter(ecg_signal[i]);
        for (int i = len - 1; i > -1; --i)
            baseline[i] = baseline_filter_.filter(baseline[i]);
        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = bandpass_filter_.filter(ecg_signal[i]);
        for (int i = len - 1; i > -1; --i)
            filt_signal[i] = bandpass_filter_.filter(ecg_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i] * filt_signal[i]);
        for (int i = len - 1; i > -1; --i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            threshold_signal[i] = threshold_filter_.filter(filt_signal[i]);
        for (int i = len - 1; i > -1; --i)
            threshold_signal[i] = threshold_filter_.filter(filt_signal[i]);

        for (unsigned int i = 0; i < len; ++i)
        {
            if (searching_for_peaks_ && (filt_signal[i] > threshold_signal[i] * threshold_ratio_) && (previous_sig_val_ > filt_signal[i]))
            {
                if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
                {
                    previous_peak_amplitude_ = previous_sig_val_;
                    samples_after_peak_count_ = 1;
                    searching_for_peaks_ = false;
                }
                else
                    previous_peak_amplitude_ *= peak_attenuation_;
            }
            else if (previous_sig_val_ < filt_signal[i])
            {
                searching_for_peaks_ = true;
                samples_after_peak_count_ = 0;
            }

            previous_sig_val_ = filt_signal[i];

            if (samples_after_peak_count_)
                ++samples_after_peak_count_;

            if (samples_after_peak_count_ == nr_slope_samples_)
            {
                samples_after_peak_count_ = 0;
                peak_signal[i] = ((marker_val_ == -1.0) ? filt_signal[i] : marker_val_);
            }
            else
                peak_signal[i] = 0;
        }
        unsigned int nr_peaks = 0;
        for (unsigned int i = nr_slope_samples_; i < len; ++i)
            if (peak_signal[i])
            {
                peak_signal[i - nr_slope_samples_ + 1] = peak_signal[i];
                peak_signal[i] = 0;
                ++nr_peaks;
            }
        const int radius = (10.0 * sampling_rate_) / 1000.0;
        for (unsigned int i = radius; i < len - radius; ++i)
            if (peak_signal[i])
            {
                unsigned int maxindx = 0, minindx = 0;
                double maxval = -2000000, minval = 2000000;
                for (int j = -radius; j < radius; ++j)
                {
                    if (maxval < ecg_signal[i + j] - baseline[i + j])
                    {
                        maxval = ecg_signal[i + j] - baseline[i + j];
                        maxindx = i + j;
                    }
                    if (minval > ecg_signal[i + j] - baseline[i + j])
                    {
                        minval = ecg_signal[i + j] - baseline[i + j];
                        minindx = i + j;
                    }
                }
                double peakval = peak_signal[i];
                peak_signal[i] = 0;
                if (maxval > -minval)
                    peak_signal[maxindx] = peakval;
                else
                    peak_signal[minindx] = peakval;
            }
        if (peak_indexes)
        {
            peak_indexes->resize(nr_peaks);
            nr_peaks = 0;
            for (unsigned int i = 0; i < len; ++i)
                if (peak_signal[i])
                    (*peak_indexes)[nr_peaks++] = i;
        }
        delete[] baseline;
    }
};
