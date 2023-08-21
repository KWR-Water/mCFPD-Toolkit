# mCFPD toolkit

The 'mCFPD toolkit' is a set of Python scripts that offers multiple approaches to performing mCFPD background leakage analyses. mCFPD (multi-DMA Comparison of Flow Pattern Distributions) is a generalized version of the CFPD method that was orinally published in 2013. 

CFPD is aimed at comparing flows into a DMA between different time periods, allowing its user to identify changes and distinguish their possible origins (e.g. leakage vs. changes in legitimate demand or water theft) .

mCFPD is aimed at comparing flows into DMAs between these DMAs, allowing its user to quantify background leakage levels. 

The toolkit contains an elaborate non-graphical script and a rudimentary GUI application. 

The mCFPD method is described in:
van Thienen, P. (2022). Direct assessment of background leakage levels for individual district metered areas (DMAs) using correspondence of demand characteristics between DMAs. Water Supply, 22(7), 6370-6388.
https://iwaponline.com/ws/article/22/7/6370/89553/Direct-assessment-of-background-leakage-levels-for

The original CFPD method is described in the following publication:
van Thienen, P. (2013). A method for quantitative discrimination in flow pattern evolution of water distribution supply areas with interpretation in terms of demand and leakage. Journal of Hydroinformatics, 15(1), 86-102.
https://iwaponline.com/jh/article/15/1/86/3218/A-method-for-quantitative-discrimination-in-flow

A version of the mCFPD script was used in the analyses presented in the 2022 paper. Development of the GUI application was started to facilitate use of the method by students. 
With shifting priorities, the GUI application was not completed to the level that was initially foreseen, but still it is fully functioning, if rudimentary. 
In my opinion, there is a lot more to be done in exploring the potential and practical applicability of the mCFPD method and developing more functionality in the toolkit.

## Installation

Download the source code and install the required dependencies. 

## Dependencies

The mCFPD Toolkit uses the following non-standard libraries:

- numpy 
- matplotlib
- scipy
- wx
- pygame

## Usage

### GUI application
Running wxclustCFPD.py launches the rudimentary GUI application. Multiple data files can be read through the File -> Load data menu (repeatedly). 
Analysis time windows can be selected either by clicking and holding on the time series plot, or by clicking "Match", which will find the best matching period between all included datasets, which can be applied by pressing "Apply".
This selected time window is marked by a cyan box. 
The user can set the comparison period length, number of bars in mCFPD diagrams, and maximum number of subclusters in the interface.
Clicking "mCFPD" will perform an mCFPD analysis for the selected datasets and settings.

The user can mark two additional time windows in each plot by using CTRL-click and SHIFT-click when selecting a time windows. These will be marked by green and blue boxes. A traditional CFPD analysis can be performed between any of the selected time windows (blue/green), also between datasets. 

### Script

Using the clusterCFPD2.py script requires editing of parameters in single() and multi() functions (for single and multiple dataset analyses, respectively), and selecting which one to call in main().

The main parameters to set are the following:
  * perlen: the length of the comparion period (in days)
  * nbar=4: the number of bars that is shown in the mCFPD diagrams, i.e. the number of scenarios between the lowest and highst plausible total background leakage level
  * maxnclust: the maximum number of coherent subclusters within the group of DMAs to consider
  
The output of the script is written in subfolders of ./output

## License

`mCFPD Toolkit` is available under a GPL3 license (https://www.gnu.org/licenses/gpl-3.0.html). 

## Development and contributions

The code was written in 2021-2022. It is not actively developed or maintained by the original developer because of a shifted work focus and lack of time. However, further cleaning up of the code, development and application are encouraged! 

## Citing

If you publish work based on `mCFPD Toolkit`, I would appreciate a citation of the following reference:

van Thienen, P. (2022). Direct assessment of background leakage levels for individual district metered areas (DMAs) using correspondence of demand characteristics between DMAs. Water Supply, 22(7), 6370-6388.
https://iwaponline.com/ws/article/22/7/6370/89553/Direct-assessment-of-background-leakage-levels-for