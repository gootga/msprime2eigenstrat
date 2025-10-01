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
   Pop1  Pop2 Pop3  Pop4  Pop5  Pop6         Pop7    Mix      Pop8       Pop9  Pop10  Pop11  Pop12          Pop13  Pop14

#Example usage: python3 msprime2eigenstrat.py
