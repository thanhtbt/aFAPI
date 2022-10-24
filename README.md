# aFAPI: Robust Subspace Tracking with Contamination Mitigation via Alpha-Divergence

We studied the problem of robust subspace tracking (RST) in contaminated environments. Leveraging the fast approximated power iteration and α-divergence, a novel robust algorithm called αFAPI was developed for tracking the underlying principal subspace of streaming data over time. αFAPI is fast and
it outperforms many RST methods while only having a low complexity linear to the data dimension. 


## DEMO

+ Run "demo_xyz.m" for synthetic experiments.

## State-of-the-art algorithms for comparison

+ **FAPI**: Badeau, Roland, Bertrand David, and Gaël Richard. "Fast approximated power iteration subspace tracking." IEEE Trans. Signal Process. 53.8 (2005): 2931-2941.

+ **LORAF**: Strobach, Peter. "The fast recursive row-Householder subspace tracking algorithm." Signal Process. 89.12 (2009): 2514-2528.

+ **GYAST**: Arjomandi-Lari, Mostafa, and Mahmood Karimi. "Generalized YAST algorithm for signal subspace tracking." Signal Process. 117 (2015): 82-95.

+ **ROBUSTA**: Linh-Trung, Nguyen, et al. "Low-complexity adaptive algorithms for robust subspace tracking." IEEE J. Sel. Topics Signal Process. 12.6 (2018): 1197-1212.

+ **RYAST**: Nguyen, Viet-Dung, Nguyen Linh Trung, and Karim Abed-Meraim. "Robust subspace tracking algorithms using fast adaptive Mahalanobis distance." Signal Process. 195 (2022): 108402.


+ **TRPAST**: A. M. Rekavandi, A.-K. Seghouane, and K. Abed-Meraim, “TRPAST: A tunable and robust projection approximation subspace tracking method,” IEEE Trans. Signal Process. (2022). 


## Results

+ Effect of p and alpha

![1 0](https://user-images.githubusercontent.com/26319211/197474765-75e53a21-1333-4c33-9595-460663db1d69.PNG)


+ Tracking in contaminated environment

![1 1](https://user-images.githubusercontent.com/26319211/197474872-10db5b21-8841-4a2f-9320-51e739596733.PNG)


+ DOA Tracking

![1 2](https://user-images.githubusercontent.com/26319211/197474904-85e9da77-542e-4d45-bc99-644e66222aba.PNG)

## Reference

This code is free and open source for research purposes. If you use this code, please acknowledge the following paper.

[1] **L.T. Thanh**, A.M. Rekavandi, S. Abd-Krim, & K. Abed-Meraim. “[*Robust Subspace Tracking with Contamination Mitigation via Alpha-Divergence*](https://)”. **TechRxiv**, 2022. [[PDF]](https://thanhtbt.github.io/files/2023_aFAPI.pdf).
