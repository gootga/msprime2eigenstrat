# msprime2eigenstrat
#Very simple script to simulate an msprime demography (based on one diploid chromosome) and convert to eigenstrat format. #Plots a phylogeny and a writes out a questionable tree as intermediate step.

#This is the actual phylogeny/demography. Pop15 (Mix) is an admixture of Pop7 (20%) and Pop8 (80%).
<img width="877" height="313" alt="image" src="https://github.com/user-attachments/assets/60d07d30-805d-48d4-991a-018248dc2f2a" />


                                                                                          ANC
                                                                   ┌───────────────────────┴─────────────────────┐
                                                                   │                                             │
                                                  ┌────────────────┴────────────────┐                            │
                                                  │                                 │                            │
                                                  │                                 │                            │
                                   ┌──────────────┴────────────┐                    │                            │
                                   │                           │                    │                            │
                       ┌───────────┴──────────┐                │                    │                            │
                       │                      │                │                    │                            │
             ┌─────────┴─────┐                │                │                    │                            │
             │               │                │                │                    │                            │
        ┌────┴────┐          │                │                │              ┌─────┴─────┐                      │
        │         │          │                │----------------│              │           │                      │
     ┌──┴──┐   ┌──┴──┐    ┌──┴──┐             │       │        │          ┌───┴───┐   ┌───┴───┐               ┌──┴──┐ 
     │     │   │     │    │     │             │       │        │          │       │   │       │               │     │ 
     P1    P2  P3    P4  P5    P6            P7      Mix      P8          P9    P10   P11    P12              P13    P14

#Example usage: python3 msprime2eigenstrat.py
