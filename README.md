# OTFS Channel Estimation

### Introduction
This project implements a **SISO-OTFS communication system** considering **fractional delay** and **fractional Doppler**, along with an **embedded single-pilot-based OTFS channel estimation** method.

---

### Delay-Doppler (DD) Domain Channel Model
The input-output relation in the DD domain is given as:

$$
y[k, l] = \sum_{k^{\prime}=0}^{N-1}\sum_{l^{\prime}=0}^{M-1} x[k^{\prime}, l^{\prime}] \, h_w[k-k^{\prime}, l-l^{\prime}] + z[k, l]
$$

- $x[k^{\prime}, l^{\prime}]$: modulation symbol at grid $[k^{\prime}, l^{\prime}]$  
- $y[k, l]$: received signal sample  
- $z[k, l]$: additive noise, e.g. $z[k, l] \sim \mathcal{CN}(0,\sigma^2)$ under rectangular window  

The original DD domain channel response is:

$$
h(\nu,\tau) = \sum_{i=1}^{P} h_i \, \delta(\nu-\nu_i)\delta(\tau-\tau_i)
$$

where $P$ is the number of channel paths, $h_i$ is the path gain, $\tau_i$ is the delay, and $\nu_i$ is the Doppler shift.

---

### Effective DD Domain Channel
The effective channel is expressed as:

$$
h_w[k-k^{\prime}, l-l^{\prime}] = \sum_{i=1}^{P} \tilde{h}_i \, w(k-k^{\prime}-k_{\nu_i}, l-l^{\prime}-l_{\tau_i})
$$

with  
- $k_{\nu_i} = \nu_i N T$  
- $l_{\tau_i} = \tau_i M \Delta_f$  
- $\tilde{h}_i = h_i e^{-j2\pi \nu_i l_{\tau_i}}$  

For ideal pulse shaping and rectangular windows, the sampling function factorizes:

$$
w(k-k^{\prime}-k_{\nu}, l-l^{\prime}-l_{\tau}) = w_{\nu}(k-k^{\prime}-k_{\nu}) \cdot w_{\tau}(l-l^{\prime}-l_{\tau})
$$

where

$$
w_\nu(\Delta k) = \frac{1}{N} e^{-j\pi (N-1)\frac{\Delta k}{N}} \frac{\sin(\pi \Delta k)}{\sin(\pi \Delta k / N)}
$$

$$
w_\tau(\Delta l) = \frac{1}{M} e^{-j\pi (M-1)\frac{\Delta l}{M}} \frac{\sin(\pi \Delta l)}{\sin(\pi \Delta l / M)}
$$

---

### Threshold-Based Channel Estimation
The proposed method estimates the **original DD domain channel response** $h(\nu,\tau)$, then reconstructs the effective channel.  
For details, see:  
P. Raviteja, K. T. Phan, and Y. Hong, *“Embedded pilot-aided channel estimation for OTFS in delay–Doppler channels,”* IEEE Trans. Vehicular Technology, 2019.

---

### Code Modules
- **Main function**  
  - `main.m`: Entry point, simulates DD domain transmission, pilot insertion, and channel estimation  

- **DD domain channel functions**  
  - `DD_Channel_Output.m`: Output matrix expression  
  - `h_omega.m`: Effective DD domain channel response  
  - `omega_nu.m`: Doppler sampling function  
  - `omega_tau.m`: Delay sampling function  

- **Channel estimation**  
  - `Insert_Pilot.m`: Pilot insertion  
  - `martix_for_CE.m`: Extract pilot info at receiver  
  - `Yu_zhi_est.m`: Threshold-based estimation  
  - `get_parameters.m`: Recover parameters from estimated response  
  - `calculate_NMSE.m`: Compute normalized mean square error (NMSE)  

---
