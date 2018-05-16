
TYPE: Unstructured
NUMBER OF BLOCKS: N/A
NUMBER OF ELEMENTS: 81,000 (approx)
TYPES OF ELEMENTS: Tet
PLUME: No
APPLICABLE DATASETS: Data8
NOTES: This model is exactly the same as the model in Grid2, except for the 
addition of several new BCs.  Each area of the model now has it's own BC so 
that symmetrix remeshing software doesn't get confused.  The two solid surface 
rings at each end of the propellant have BCs called SlideInXHeadEnd and 
SlideInXAftEnd.  The solid surface patches at the head end and nozzle of the 
rocket now have separate BCs.  And the source area and farfield plume BCs are
left the same.  

8/2/06: Added new BC called ErodingNozzle for modeling of 2-part Rocburn for
erosion of throat. 

8/22/06: Removed plume.  Otherwise the same as Grid5. 

02/27/08: Significantly coarsened, also meshed the head, aft, and center
          bore separately in an attempt to do more significant burn back
          without remeshing.