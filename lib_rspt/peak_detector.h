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
    const double previous_peak_reference_attenuation_ = 15;
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
    inline double detect(double new_sample)
    {
        if (!sample_indx_++)
            bandpass_filter_.init_history_values(new_sample, sampling_rate_);

        double sig_val = bandpass_filter_.filter(new_sample);
        sig_val = integrative_filter_.filter(sig_val * sig_val);
        double threshold = threshold_filter_.filter(sig_val);

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
