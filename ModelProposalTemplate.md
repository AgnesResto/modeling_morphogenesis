# Model Proposal for _[Project Name]_

_Your Name_

* Course ID: CMPLXSYS 530,
* Course Title: Computer Modeling of Complex Systems
* Term: Winter, 2019



&nbsp; 

### Goal 
*****
 
_The goal of this human embryonic stem cell (hESC) morphogenesis ABM will be to understand hESC self-organization and differentiation dynamics in different in vitro environments._

&nbsp;  
### Justification
****
_Laboratory experiments with hESCs often require vast amounts of resources in terms of both time and money. For this reason, I propose the development of an in silico model to test hypothesis regarding morphogenetic events of hESCs. Because of the discrete nature of cells, an agent-based model (ABM) will be used to model these events._

&nbsp; 
### Main Micro-level Processes and Macro-level Dynamics of Interest
****

_The model will seek to understand two key processes of hESC morphogenesis: (1) cyst formation, and (2) hESC differentiation. For the cyst formation, the model will be looking at how a cell acts according to its local neighborhood and how this can bring about self-organization in the system that leads to the formation of spherical polarized cysts. With regards to differentiation, the model will be looking at how contact with matrix in the environment can lead to a cell’s differentiation. Additionally, the model will explore the possible inductive effects of differentiated cells (i.e. could differentiated cells be secreting signals that make other cells differentiate). The model will also look at possible inhibition of differentiation by pluripotent cells._

&nbsp; 


## Model Outline
****
&nbsp; 
### 1) Environment
_Description of the environment in your model. Things to specify *if they apply*:_

* _Boundary conditions (e.g. wrapping, infinite, etc.)_
* _Dimensionality (e.g. 1D, 2D, etc.)_
* _List of environment-owned variables (e.g. resources, states, roughness)_
* _List of environment-owned methods/procedures (e.g. resource production, state change, etc.)_

_The model will be constructed in Netlogo in a 2D hexagonal grid with no wrapping._
Environment-owned variables:
1.	Amount of diffused inhibition chemical from pluripotent stem cells
2.	Amount of diffused induction chemical from differentiated cells
3.	Amount of matrix in the environment
Environment-owned procedures:
1.	Diffusion of chemical in the environment


```python
# Include first pass of the code you are thinking of using to construct your environment
# This may be a set of "patches-own" variables and a command in the "setup" procedure, a list, an array, or Class constructor
# Feel free to include any patch methods/procedures you have. Filling in with pseudocode is ok! 
# NOTE: If using Netlogo, remove "python" from the markdown at the top of this section to get a generic code block
patches-own [
  patch-neighbors-6
  inhibitor
  induction-chemical

]

to setup
  clear-all
  if culture-condition = "embedded"
  [make-embedded-culture]
  if culture-condition = "clustered"
  [make-clusters]
  make-clusters
  ;;make-embedded-culture
  ask cells [set can-divide? True]
  ask cells [set counted? False]
  ask cells [set num-divisions 0]
  set ticks-matrix-differentiation 0
  set lumen-num 0
  set num-cysts 1
  set steady-state? False
  reset-ticks
end


to make-embedded-culture
  set-default-shape matrix "box"
  ask n-of num-cells patches [
    sprout 1
    [set breed cells
      set ss? False
      set shape atom-shape
      set color red
      if pxcor mod 2 = 0 [set ycor ycor - 0.5 ]
    ]
  ]
  ask patches [
    if not any? turtles-here [
    if overlay-matrix-percentage > random 100 [
    sprout 1 [
    set breed matrix
    set color white
  if pxcor mod 2 = 0 [set ycor ycor - 0.5 ]
      ]
      ]
    ]
  ]
  ask cells [
     ;; define neighborhood of patches
     ifelse pxcor mod 2 = 0 [
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]
end

to make-clusters
  set-default-shape matrix "box"
  ask patches [
    if not any? turtles-here [
      if overlay-matrix-percentage > random 100 [
        sprout 1 [
          set breed matrix
          set color white
          if pxcor mod 2 = 0 [set ycor ycor - 0.5 ]
        ]
      ]
    ]
  ]
  while [count cells < num-cells]
  [
    ask n-of 3 matrix [
      if (ycor > min-pycor + 5) and (ycor < max-pycor - 5) and (xcor > min-pxcor + 5) and (xcor < max-pxcor - 5)
      [
        hatch-cells 1
        [
          set breed cells
          set color red
          set shape atom-shape
          set ss? False
        ]
        ifelse cluster-size > 1 or cluster-size < 1
        [
          ask matrix in-radius cluster-size
          [
            hatch-cells 1
            [
              set breed cells
              set color red
              set shape atom-shape
              set ss? False
            ]
          ]
        ]
        [
          if cluster-size = 1 [
            ask n-of 2 matrix in-radius 1.5
            [
              hatch-cells 1
              [
                set breed cells
                set color red
                set shape atom-shape
                set ss? False
              ]
            ]
          ]
        ]
      ]
    ]
  ]
  ask cells [
    if any? matrix-here [ask matrix-here [die]]
    if count cells-here > 1 [ask one-of cells-here [die]]
    ;; define neighborhood of patches
    ifelse pxcor mod 2 = 0 [
      set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
      set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]
end


```

&nbsp; 

### 2) Agents
 
 _Description of the "agents" in the system. Things to specify *if they apply*:_
 
* _List of agent-owned variables (e.g. age, heading, ID, etc.)_
* _List of agent-owned methods/procedures (e.g. move, consume, reproduce, die, etc.)_

Agent	Variables	Procedures
Cells	
neighbors-6: stores information about the type of neighboring agents and their locations
cycles-matrix-contact: how many ticks has a cell been in contact with matrix
cycles-diff-contace: how many ticks has a cell been in contact with a differentiated cell
still_count: how many cycles has a cell stood still
ss?: has the cell reached steady state?
can-divide?: can the cell divide?
group-id: identifies to which cyst the cell belongs to
num-divisions: keeps track of how many times a cell has divided
A cell can move, die, divide, consume matrix, make links to other cells, and differentiate
Matrix	
Agents that are a part of the environment	Matrix can be consumed by cells and cause cells to differentiate
Lumen	
These agents are used to identify the inner cavity of cysts	N/A
Walker	
This agent is created once the cysts form and is used to give a unique identifier to the lumen and cells of each cyst	Move, assign id to lumen from a cyst
Cell-links	
These are links formed between two cells that indicate the start of polarization.	N/A

```python
# Include first pass of the code you are thinking of using to construct your agents
# This may be a set of "turtle-own" variables and a command in the "setup" procedure, a list, an array, or Class constructor
# Feel free to include any agent methods/procedures you have so far. Filling in with pseudocode is ok! 
# NOTE: If using Netlogo, remove "python" from the markdown at the top of this section to get a generic code block
breed [ cells cell ]
breed [ matrix a-matrix ]
breed [ glass a-glass ]
breed [ lumen a-lumen ]
breed [ walker a-walker ]
undirected-link-breed [ cell-links cell-link ]
undirected-link-breed [ neighbor-links neighbor-link ]
undirected-link-breed [ lumen-links lumen-link ]
;breed [ differentiated diff ]

cells-own [
  neighbors-6
  cycles-matrix-contact ;;counts number of cycles in contact with >= x amount of     matrix neighbors
  cycles-diff-contact ;;counts number of cycles in contact with >= x amount of differentiated neighbors
  still_count ;; counts number of cycles a cell has kept still
  ss? ;;boolean stating if cell has reached steady state
  counted?
  can-divide?
  group-id
  num-divisions
]

lumen-own
[lumen-neighbors-6
 lumen-id
 group-id
 active?
 unvisited?
 gateway?
 gateway-visited?
]

walker-own [
 id-in-proximity?
 walker-id
 up-moves
 down-moves
 group-id
]

matrix-own [
  matrix-neighbors-6
]

# The rest of the code can be found under morphogenesis. There are three working versions: Morphogenesis_3Denv_3Dov, Morphogenesis_3Denv_3Dov-v2, and Morphogenesis_3Denv_3Dov-v3



```

&nbsp; 

### 3) Action and Interaction 
 
**_Interaction Topology_**

_The interaction topology is most similar to a CA neighborhood. The cells will only take actions and move within their local neighborhood. _
 
**_Action Sequence_**

_What does an agent, cell, etc. do on a given turn? Provide a step-by-step description of what happens on a given turn for each part of your model_

1.	The cells check their variables with regards to differentiation (i.e. do they differentiate or stay pluripotent)
2.	The cells check their local neighborhood and decide to move, die, consume matrix, divide, or differentiate.
3.	The cells update their local environment.
4.	Once the system has reached stead state (i.e. the cells have established satisfaction with their state), the lumen of each cyst is identified with a walker that moves through the environment giving unique id to all the lumen of each cyst. After the lumen have their IDs they pass it on to the cells of the cyst.


&nbsp; 
### 4) Model Parameters and Initialization

Global variables:
1.	ticks-matrix-differentiation- time in which matrix can cause differentiation
2.	steady-state?- has the system reached steady state
3.	num-cysts- total number of cysts in the system
4.	grid-x-pos- x coordinate the walker is on
5.	grid-y-pos- y coordinate the walker is on

The model setup includes creating an environment that is filled with a percentage of matrix (chosen by the model user) and positioning a number of cells (chosen by the model user) in either a random or clustered manner. The model user has to establish the following parameters: amount of cells, cluster size, maximum number of divisions for each cell and whether there will be diffusion in the environment or differentiation will take place through contact. If there is diffusion, the model user must also establish how much chemical (both from induction and inhibition) the cells secrete, and what are the thresholds for induction and inhibition. If there is no diffusion, the model user must establish how many matrix neighbors will lead to differentiation, how many differentiated neighbors will lead to differentiation, and the respective number of cycles through which the contact must be maintained before differentiation. 
Schedule during each tick:

1.	Chemicals diffuse in the environment
2.	Chemicals degrade in the environment
3.	Cells check their state with regards to differentiation
4.	Cells release either inhibitor or inductive signaling into the environment
5.	Cells check their neighbors and take an action (die, move, divide, consume matrix)
6.	Cells update their neighborhood
7.	If the system has reached steady state:
a.	Patches establish their neighborhood
b.	Lumen is created in the cavities of cysts
c.	Walkers are created and give unique IDs to lumen
d.	Lumen pass on their IDs to the cells in the cyst
8.	If all the cells in the system have group-ids the model stops

&nbsp; 

### 5) Assessment and Outcome Measures

_What quantitative metrics and/or qualitative features will you use to assess your model outcomes?_

 Model calibration will be performed for: (1) Cyst formation, and (2) Cell differentiation. In terms of cyst formation, I will be looking at the number of cells in each cyst and the cyst shape. In terms of cell differentiation, I will be looking for model parameters that give pluripotent cysts at high plating densities coupled with differentiated cysts at low plating densities. One particular even of interest is the formation of asymmetric cysts. This kind of event is rare in laboratory experiments, and the in silico model will be used to test hypothesis with regards to how it occurs.

&nbsp; 

### 6) Parameter Sweep

_What parameters are you most interested in sweeping through? What value ranges do you expect to look at for your analysis?_

The parameters I will be varying initially are:
Differentiation through contact
number of matrix contacts required for differentiation	1-5
cycles of contact with matrix for differentiation	1-10
number of differentiated neighbors required for differentiation	1-5
cycles of contact with differentiated cells for differentiation	1-10
Diffusion in the system
amount of inhibitor released by pluripotent cells	undetermined
amount of inductive chemical released by differentiated cells	undetermined
Induction chemical threshold for differentiation	undetermined
maximum number of divisions per cell	1-3
*With regards to examining diffusion, I’ll be looking at the relationship between the parameters rather than actual quantitative values. 
