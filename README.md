# alankylam.github.io
# Data Science Portfolio by Alan Lam

This portfolio is a compilation of projects which I created for corporate data analysis and undergraduate courses. 

## Corporate Drilling Data Project 

This project tackled the data integrity issue that caused major discrepancies between theoretical and actual well data. 

There were major differences in the field captured and software calculated bottom-holes (BH), the deepest known point of a well, within the corporate repository of over 19,000 wells. 

[ Insert Image of Well] 
[Insert Image of discrepenacy] 
### Analysis
First, we must understand how drilling data is recorded and visualized. 

Drilling data is collected through measurement-while-drilling (MWD) tools during the drilling process that utilizes electrical impulses to transmit a combination of data. The survey data collected typically consists of measured depth, inclination, and azimuth--or MDINCAZM for short. 
[Insert drilling gif]
The MDINCAZM data is then visualized using the mathematical formula called the minimum curvature method displayed below.
[insert MD picture]
With this in mind, my analysis first focused on the different parts of the well to understand which part of had the largest impact.
[insert analysis 1, different parts of the well]
The conclusion from this analysis found that the vertical length was the major culprit. The vertical section of the data set was left out for operational cutbacks in order to decrease the cost of the MWD tool use. But, in return, this caused data integrity issues for the office. 
[insert missing section grapgic] [with missing section data graphic]
I conducted further analysis on the missing vertical lengths to understand on what magnitude this impacted the BH discrepancy. 
[insert analysis 2, analysis of missing length]
The conclusion from this was that a larger missing vertical data set for a well did not directly lead to a larger discrepancy, but rather, a larger error potential.  
To take my analysis further, I was then able to create a heatmap that highlighted BH errors alongside well parameters to further understand this discrepancy. 
[insert heatmap photos]
### Solution
Next, I began testing for different solutions to minimize the BH difference. 
I began interpolating the missing well data to fill in the data gaps for various wells. After defining interpolation rules, I was able to minimize the BH discrepancy by ~70% (e.g. 50ft BH error to 15ft error).
[insert fixed data sheet] [insert data results]
I was then able to create a Python script to streamline this process for the corporate repository. 
[corporate repository before/after]
