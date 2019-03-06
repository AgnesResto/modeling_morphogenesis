breed [ cells cell ]
breed [ matrix matri ]
breed [ glass glas ]
breed [ lumen a-lumen ]
breed [ walker a-walker ]
undirected-link-breed [ cell-links cell-link ]
undirected-link-breed [ neighbor-links neighbor-link ]
undirected-link-breed [ lumen-links lumen-link ]
;breed [ differentiated diff ]




cells-own [
  neighbors-6
  cycles-matrix-contact ;;counts number of cycles in contact with >= x amount of matrix neighbors
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

patches-own [
  patch-neighbors-6

]

;differentiated-own [
;  neighbors-6
;]

globals [
  logtime                          ;; log of time
  colors                           ;; used both to color turtles, and for histogram
  xmax                             ;; max x size
  ymax                             ;; max y size
  ticks-matrix-differentiation     ;; time in which matrix can cause differentiation
  steady-state?
  lumen-num
  num-cysts
  grid-x-pos
  grid-y-pos
  any-lumen-neighbors?
]

to setup
  random-seed 89942
  ;;random-seed 22222
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


to go
  ;;update counter for number of cycles in contact with matrix
  ifelse steady-state? = False
  [

  ask cells [
    if diff-matrix? = True and diff-induction? = False [
    ifelse count matrix-on neighbors-6 >= num-matrix-diff
      [set cycles-matrix-contact cycles-matrix-contact + 1]
      [ifelse cycles-matrix-contact = 0
        [set cycles-matrix-contact 0]
        [set cycles-matrix-contact cycles-matrix-contact - 1]
    ]
    if ticks-matrix-differentiation < max-ticks-differentiation [
    if color = red [
    if cycles-matrix-contact >= cycles-diff-matrix [
      set color green
      set ticks-matrix-differentiation ticks-matrix-differentiation + 1]
    ]
  ]
  ]
    if diff-matrix? = True and diff-induction? = True [
      ifelse count matrix-on neighbors-6 >= num-matrix-diff
        [set cycles-matrix-contact cycles-matrix-contact + 1]
        [ifelse cycles-matrix-contact = 0
          [set cycles-matrix-contact 0]
          [set cycles-matrix-contact cycles-matrix-contact - 1]
        ]
        if color = red [
          ifelse cycles-matrix-contact >= cycles-diff-matrix
          [set color green]
          [set ticks-matrix-differentiation ticks-matrix-differentiation + 1]
        ]
      ifelse count cells-on neighbors-6 with [color = green] >= num-diff-ind and color = red
      [set cycles-diff-contact cycles-diff-contact + 1
      if cycles-diff-contact >= cycles-diff-ind and (count cells in-radius 1.5 with [color = red]) < undiff-num-inhibition
        [set color green]]
      [if color = red and cycles-diff-contact > 0
        [set cycles-diff-contact cycles-diff-contact - 1]
      ]
    ]
  ]


    ifelse rule-set = "Original"
    [
      ask cells [assess-state-even]
      ;;update neighbors after each iteration
      ask cells [
        if count cells-here > 1 [ ask one-of cells-here [die]]]
      ask cells [
        ;; define neighborhood of patches
        ifelse pxcor mod 2 = 0 [
          set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
        ][
          set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
      ]
      ask cells [assess-state-odd]
      ask cells [
        if count cells-here > 1 [ ask one-of cells-here [die]]]
      ask cells [
        ;; define neighborhood of patches
        ifelse pxcor mod 2 = 0 [
          set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
        ][
          set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
      ]
    ]
    [
      if rule-set = "New"
      [
        ask cells [
        new-rules
        ifelse pxcor mod 2 = 0 [
          set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
        ][
          set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
        ]
      ]
;      ask cells [
;        ifelse pxcor mod 2 = 0 [
;          set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
;        ][
;          set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
;      ]
    ]

    ;;create-cell-links-with cells-on neighbors-6 with [not link-neighbor? myself] [set color white]
  if (all? cells with [ycor < (max-pycor - 2) and ycor > (min-pycor + 2) and xcor > (min-pxcor + 2) and xcor < (max-pycor - 2)] [ss? = True]) [

;    ask cells [
;      ifelse pxcor mod 2 = 0 [
;        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
;    ][
;        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
;    create-cell-links-with cells-on neighbors-6 with [not link-neighbor? myself] [set color white]
;  ]
      ask patches [
     ;; define neighborhood of patches
     ifelse pxcor mod 2 = 0 [
        set patch-neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set patch-neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]
    ask matrix [
     ;; define neighborhood of patches
     ifelse pxcor mod 2 = 0 [
        set matrix-neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set matrix-neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]

;    find-cysts
;    find-lumen
    set steady-state? True
    ;ask patches [percolate3]
;   count-cysts
;    stop
  ]
]
  [if steady-state? = True [
    find-cysts
    find-lumen
    make-cyst-id
    make-cell-id
    propagate-id
    if all? lumen [group-id > 0]
    [stop]
    ]
  ]
  tick
end

to refresh-neighbors
  ifelse pxcor mod 2 = 0 [
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to new-rules
;  ifelse pxcor mod 2 = 0 [
;    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
;  ][
;    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  if (pycor > min-pycor) and (pycor < max-pycor) and (pxcor > min-pxcor) and (pxcor < max-pxcor)
  [
    ifelse count cells-here = 1
    [
      ifelse not any? turtles-on neighbors-6
      [
        die
      ]
      [
        ifelse not any? cells-on neighbors-6
        [
          ifelse count matrix-on neighbors-6 = 6
          [
            hatch-matrix 1 [set color white]
            die
            ;; this is assuming that previous connection with a cell would have led to empty patch nearby
          ]
          [
            ifelse count matrix-on neighbors-6 <= 5
            [
              ifelse not any? cell-links
              [
                die
              ]
              [
                ifelse can-divide? = True
                [
                  set num-divisions num-divisions + 1
                  if num-divisions >= max-divisions [set can-divide? False]
                  hatch-cells 1
                  [
                    create-cell-link-with myself [set color blue]
                    ifelse pxcor mod 2 = 0
                    [rule3-even]
                    [rule3-odd]
                    ifelse pxcor mod 2 = 0 [
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                    ][
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
                  ]
                ]
                [
                  rule7 ;; move to matrix close to other cells
                ]
              ]
            ]
            [
              ifelse count matrix-on neighbors-6 = 4
              [
                ;; add if necessary
              ]
              [
                ;; add if count matrix-on neighbors-6 = 3 if necessary
              ]
            ]
          ]
        ]
        [
          ifelse count cells-on neighbors-6 = 6
          [ ;; maybe having all cell neighbors is ok for a certain amount of time
            ;            ifelse any? patches in-radius 2.5 with [not any? cells-here]
            ;            [
            ;              move-to one-of patches in-radius 2.5 with [not any? cells-here and count cells in-radius 1.5 >= 1]
            ;              if pxcor mod 2 = 0 [ set ycor ycor - 0.5 ]
            ;            ]
            ;            [
            ;              die
            ;            ]
            ifelse not any? link-neighbors
            [
              push-cell-to-matrix
              ask cells[
                ifelse pxcor mod 2 = 0 [
                  set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                ][
                  set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
              ]
            ]
            [
              ;; what happens when the cell is surrounded by 6 cells but already established some links before?
            ]
          ]
          [
            ifelse count cells-on neighbors-6 = 1
            [
              ifelse count matrix-on neighbors-6 = 5
              [
                ifelse (cell-link-neighbor? (one-of cells-on neighbors-6) = False) and (not any? (cells-on neighbors-6) with [count cells-on neighbors-6 > 1]) ;; and (not any? ((cells-on neighbors-6) with [cell-link-neighbor? (one-of cells-on neighbors-6) = True]))
                  [
                    create-cell-links-with cells-on neighbors-6 [set color blue]
                    ifelse pxcor mod 2 = 0
                    [rule2-even]
                    [rule2-odd]
                ]
                [
                  ;; what happens when your neighbor is already linked to another cell?
                ]
                ifelse pxcor mod 2 = 0 [
                  set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                ][
                  set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
              ]
              [
                ifelse count matrix-on neighbors-6 = 4
                [
                  ifelse can-divide? = True
                  [
                    set num-divisions num-divisions + 1
                    if num-divisions >= max-divisions [set can-divide? False]
                    hatch-cells 1
                    [
                      create-cell-link-with myself [set color blue]
                      ifelse pxcor mod 2 = 0
                      [rule3-even]
                      [rule3-odd]
                      ifelse pxcor mod 2 = 0 [
                        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                      ][
                        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
                    ]
                  ]
                  [
                    ;; what happens when it can't divide? maybe move to one of matrix next to cell
                    hatch-matrix 1 [set color white]
                    ifelse pxcor mod 2 = 0
                    [rule10-even]
                    [rule10-odd]
                    ifelse pxcor mod 2 = 0 [
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                    ][
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
                  ]
                  ifelse pxcor mod 2 = 0 [
                    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                  ][
                    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
                ]
                [
                  ifelse count matrix-on neighbors-6 = 0
                  [
                    ;                    ifelse pxcor mod 2 = 0
                    ;                    [rule6-even]
                    ;                    [rule6-odd]
                    rule9
                  ]
                  [
                    ifelse count matrix-on neighbors-6 = 1
                    [
                      if any? matrix-on neighbors-6 with [count cells in-radius 1.5 > 1]
                      [
                        move-to one-of matrix-on neighbors-6
                        ask matrix-here [die]
                      ]
                    ]
                    [
                      ;;ifelse count matrix-on neighbors-6 = 2, 3
                    ]
                  ]
                ]
              ]
            ]
            [
              ifelse count cells-on neighbors-6 <= 3 and count matrix-on neighbors-6 >= 2 and (count cells-on neighbors-6 + count matrix-on neighbors-6) < 6
              [
                ifelse can-divide? = True
                [
                  set num-divisions num-divisions + 1
                  if num-divisions >= max-divisions [set can-divide? False]
                  hatch-cells 1
                  [
                    create-cell-link-with myself [set color blue]
                    ifelse pxcor mod 2 = 0
                    [rule5-even]
                    [rule5-odd]
                    ;; what happens when there is no matrix next to a cell? Should they move next to the cell?
                    ifelse pxcor mod 2 = 0 [
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                    ][
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
                  ]
                ]
                [
                  ;; This causes broken cysts!!
                  ;; what happens when they can no longer divide? Follow rule 5 on their own?
                  ;                    ifelse pxcor mod 2 = 0
                  ;                    [rule5-even]
                  ;                    [rule5-odd]
                  ;                    ifelse pxcor mod 2 = 0 [
                  ;                      set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                  ;                    ][
                  ;                      set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
                ]

              ]
              [
                ifelse (count cells-on neighbors-6 >= 3) and count matrix-on neighbors-6 >= 1 and (count matrix-on neighbors-6 + count cells-on neighbors-6) = 6
                [
                  let my-cell-neighbors (count cells-on neighbors-6)
                  ifelse not any? link-neighbors and not any? (cells-on neighbors-6) with [count (cells-on neighbors-6) > my-cell-neighbors]
                  [
                    create-cell-links-with (cells-on neighbors-6) with [not link-neighbor? myself] [set color blue]
                    move-to one-of matrix-on neighbors-6
                    ask matrix-here [die]
                  ]
                  [
                    ;; do something else if you have links already
                  ]
                ]
                [
                  ifelse count cells-on neighbors-6 <= 4 and count matrix-on neighbors-6 = 1
                  [
                    move-to one-of matrix-on neighbors-6
                    ask matrix-here [die]
                    ifelse pxcor mod 2 = 0 [
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                    ][
                      set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]

                  ]
                  [
                    ifelse (count cells-on neighbors-6 = 2) and (count matrix-on neighbors-6 = 4)
                    [
                      ifelse not any? link-neighbors
                      [
                        create-cell-links-with (cells-on neighbors-6) with [not link-neighbor? myself] [set color blue]
                        ifelse pxcor mod 2 = 0
                        [rule2-even]
                        [rule2-odd]
                        ifelse pxcor mod 2 = 0 [
                          set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
                        ][
                          set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
                        ask cells in-radius 2 [refresh-neighbors]
                      ]
                      [
                        ask cells in-radius 1.5 [refresh-neighbors]
                        push-cell
                        ;; if its already linked do something different
                      ]
                    ]
                    [
                      ifelse count cells-on neighbors-6 = 2 and count matrix-on neighbors-6 = 0
                      [
                        rule9
                      ]
                      [
                        ;; if count cells-on ...
                      ]
                    ]
                  ]
                ]
              ]
            ]
          ]
        ]
      ]
    ]
    [
      ask one-of cells-here [die]
    ]
  ]
end

to push-cell-to-matrix
  let movable-cells (cells-on neighbors-6) with [any? matrix in-radius 1.5]
  ifelse count movable-cells >= 1
  [
    move-to one-of movable-cells
    ifelse pxcor mod 2 = 0 [
      set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
      set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
    ask one-of cells-here
    [
      move-to one-of matrix-on neighbors-6
      ask matrix-here [die]
    ]
  ]
  [
    ;; what happens when there are no movable cells in the neighbors?
  ]
  ifelse pxcor mod 2 = 0 [
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to push-cell
  let movable-cells (cells-on neighbors-6) with [any? matrix in-radius 1.5 with [count turtles in-radius 1.5 < 7]]
  ifelse count movable-cells >= 1
  [
    hatch-matrix 1 [set color white]
    move-to one-of movable-cells
    ask one-of cells-here
    [
      ifelse pxcor mod 2 = 0
      [rule3-even]
      [rule3-odd]
    ]
  ]
  [
    ifelse can-divide? = True
    [set num-divisions num-divisions + 1
      if num-divisions >= max-divisions [set can-divide? False]
      hatch-cells 1
      [
        create-cell-link-with myself [set color blue]
        ifelse pxcor mod 2 = 0
        [rule5-even]
        [rule5-odd]
        ifelse pxcor mod 2 = 0 [
          set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
        ][
          set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
      ]
    ]
    [
      ;; what happens when there's 4 matrix + 2 cells around and the cell can no longer divide?
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule3-even
  ifelse not any? turtles-on patch-at 0 1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at -1 0)
    [
      ifelse any? matrix-on patch-at 1 0
      [
        move-to patch-at 1 0
        ask matrix-here [die]
      ]
      [
        if any? matrix-on patch-at -1 0
        [
          move-to patch-at -1 0
          ask matrix-here [die]
        ]
      ]
    ]
    [
      ifelse not any? turtles-on patch-at 1 0 and (any? matrix-on patch-at 0 1 or any? matrix-on patch-at 1 -1)
        [
          ifelse any? matrix-on patch-at 0 1
          [
            move-to patch-at 0 1
            set ycor ycor - 0.5
            ask matrix-here [die]
          ]
          [
            if any? matrix-on patch-at 1 -1
            [
              move-to patch-at 1 -1
              ask matrix-here [die]
            ]
          ]
        ]
        [
          ifelse not any? turtles-on patch-at 1 -1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at 0 -1)
            [
              ifelse any? matrix-on patch-at 1 0
              [
                move-to patch-at 1 0
                ask matrix-here [die]
              ]
              [
                if any? matrix-on patch-at 0 -1
                [
                  move-to patch-at 0 -1
                  set ycor ycor - 0.5
                  ask matrix-here [die]
                ]
              ]
            ]
            [
              ifelse not any? turtles-on patch-at 0 -1 and (any? matrix-on patch-at 1 -1 or any? matrix-on patch-at -1 -1)
                [
                  ifelse any? matrix-on patch-at 1 -1
                  [
                    move-to patch-at 1 -1
                    ask matrix-here [die]
                  ]
                  [
                    if any? matrix-on patch-at -1 -1
                    [
                      move-to patch-at -1 -1
                      ask matrix-here [die]
                    ]
                  ]
                ]
                [
                  ifelse not any? turtles-on patch-at -1 -1 and (any? matrix-on patch-at 0 -1 or any? matrix-on patch-at -1 0)
                    [
                      ifelse any? matrix-on patch-at 0 -1
                      [
                        move-to patch-at 0 -1
                        set ycor ycor - 0.5
                        ask matrix-here [die]
                      ]
                      [
                        if any? matrix-on patch-at -1 0
                        [
                          move-to patch-at -1 0
                          ask matrix-here [die]
                        ]
                      ]
                    ]
                    [
                      if not any? turtles-on patch-at -1 0 and (any? matrix-on patch-at -1 -1 or any? matrix-on patch-at 0 1)
                        [
                          ifelse any? matrix-on patch-at -1 -1
                          [
                            move-to patch-at -1 -1
                            ask matrix-here [die]
                          ]
                          [
                            if any? matrix-on patch-at 0 1
                            [
                              move-to patch-at 0 1
                              set ycor ycor - 0.5
                              ask matrix-here [die]
                            ]
                          ]
                        ]
                    ]
                ]
            ]
        ]
    ]
end

to rule3-odd ;; cell neighbors <= 1 and some empty spaces - move to matrix adjacent to empty space
  ifelse not any? turtles-on patch-at 0 1 and (any? matrix-on patch-at -1 1 or any? matrix-on patch-at 1 1)
    [
      ifelse any? matrix-on patch-at -1 1
      [
        move-to patch-at -1 1
        set ycor ycor - 0.5
        ask matrix-here [die]
      ]
      [
        if any? matrix-on patch-at 1 1
        [
          move-to patch-at 1 1
          set ycor ycor - 0.5
          ask matrix-here [die]
        ]
      ]
    ]
    [
      ifelse not any? turtles-on patch-at 1 1 and (any? matrix-on patch-at 0 1 or any? matrix-on patch-at 1 0)
        [
          ifelse any? matrix-on patch-at 0 1
          [
            move-to patch-at 0 1
            ask matrix-here [die]
          ]
          [
            if any? matrix-on patch-at 1 0
            [
              move-to patch-at 1 0
              set ycor ycor - 0.5
              ask matrix-here [die]
            ]
          ]
        ]
        [
          ifelse not any? turtles-on patch-at 1 0 and (any? matrix-on patch-at 1 1 or any? matrix-on patch-at 0 -1)
            [
              ifelse any? matrix-on patch-at 1 1
              [
                move-to patch-at 1 1
                set ycor ycor - 0.5
                ask matrix-here [die]
              ]
              [
                if any? matrix-on patch-at 0 -1
                [
                  move-to patch-at 0 -1
                  ask matrix-here [die]
                ]
              ]
            ]
            [
              ifelse not any? turtles-on patch-at 0 -1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at -1 0)
                [
                  ifelse any? matrix-on patch-at 1 0
                  [
                    move-to patch-at 1 0
                    set ycor ycor - 0.5
                    ask matrix-here [die]
                  ]
                  [
                    if any? matrix-on patch-at -1 0
                    [
                      move-to patch-at -1 0
                      set ycor ycor - 0.5
                      ask matrix-here [die]
                    ]
                  ]
                ]
                [
                  ifelse not any? turtles-on patch-at -1 0 and (any? matrix-on patch-at -1 1 or any? matrix-on patch-at 0 -1)
                    [
                      ifelse any? matrix-on patch-at -1 1
                      [
                        move-to patch-at -1 1
                        set ycor ycor - 0.5
                        ask matrix-here [die]
                      ]
                      [
                        if any? matrix-on patch-at 0 -1
                        [
                          move-to patch-at 0 -1
                          ask matrix-here [die]
                        ]
                      ]
                    ]
                    [
                      if not any? turtles-on patch-at -1 1 and (any? matrix-on patch-at -1 0 or any? matrix-on patch-at 0 1)
                        [
                          ifelse any? matrix-on patch-at -1 0
                          [
                            move-to patch-at -1 0
                            set ycor ycor - 0.5
                            ask matrix-here [die]
                          ]
                          [
                            if any? matrix-on patch-at 0 1
                            [
                              move-to patch-at 0 1
                              ask matrix-here [die]
                            ]
                          ]
                        ]
                    ]
                ]
            ]
        ]
    ]
end

to rule2-even
  ifelse any? cells-on patch-at 0 1
  [
    move-to patch-at 0 -1
    set ycor ycor - 0.5
    ask matrix-here [die]
  ]
  [ifelse any? cells-on patch-at 1 0
    [
      move-to patch-at -1 -1
      ask matrix-here [die]]
    [ifelse any? cells-on patch-at 1 -1
      [
        move-to patch-at -1 0
        ask matrix-here [die]
      ]
      [ifelse any? cells-on patch-at 0 -1
        [
          move-to patch-at 0 1
          set ycor ycor - 0.5
          ask matrix-here [die]
        ]
        [ifelse any? cells-on patch-at -1 -1
          [
            move-to patch-at 1 0
            ask matrix-here [die]
          ]
          [
            move-to patch-at 1 -1
            ask matrix-here [die]
          ]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    ;; update neighborhood for cell
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule2-odd ;; only one neighbor with no other neighbors - forma a link and move opposite leaving a lumen
  ifelse any? cells-on patch-at 0 1
  [
    move-to patch-at 0 -1
    ask matrix-here [die]
  ]
  [ifelse any? cells-on patch-at 1 1
    [
      move-to patch-at -1 0
      set ycor ycor - 0.5
      ask matrix-here [die]
    ]
    [ifelse any? cells-on patch-at 1 0
      [
        move-to patch-at -1 1
        set ycor ycor - 0.5
        ask matrix-here [die]
      ]
      [ifelse any? cells-on patch-at 0 -1
        [
          move-to patch-at 0 1
          ask matrix-here [die]
        ]
        [ifelse any? cells-on patch-at -1 0
          [
            move-to patch-at 1 1
            set ycor ycor - 0.5
            ask matrix-here [die]
          ]
          [
            move-to patch-at 1 0
            set ycor ycor - 0.5
            ask matrix-here [die]
          ]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    ;; update neighborhood for cell
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule5-even
  ifelse (any? matrix-on patch-at 0 1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at -1 0))
  [
    move-to patch-at 0 1
    set ycor ycor - 0.5
    ask matrix-here [die]
  ]
  [ifelse (any? matrix-on patch-at 1 0 and (any? cells-on patch-at 0 1 or any? cells-on patch-at 1 -1))
    [
      move-to patch-at 1 0
      ask matrix-here [die]
    ]
    [ifelse (any? matrix-on patch-at 1 -1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at 0 -1))
      [
        move-to patch-at 1 -1
        ask matrix-here [die]
      ]
      [ifelse (any? matrix-on patch-at 0 -1 and (any? cells-on patch-at 1 -1 or any? cells-on patch-at -1 -1))
        [
          move-to patch-at 0 -1
          set ycor ycor - 0.5
          ask matrix-here [die]
        ]
        [ifelse (any? matrix-on patch-at -1 -1 and (any? cells-on patch-at 0 -1 or any? cells-on patch-at -1 0))
          [
            move-to patch-at -1 -1
            ask matrix-here [die]
          ]
          [if (any? matrix-on patch-at -1 0 and (any? cells-on patch-at -1 -1 or any? cells-on patch-at 0 1))
            [
              move-to patch-at -1 0
              ask matrix-here [die]
            ]
          ]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    ;; update neighborhood for cell
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule5-odd ;; if there are two neighbor cells and three matrix neighbors move to matrix next to cell
  ifelse (any? matrix-on patch-at 0 1 and (any? cells-on patch-at 1 1 or any? cells-on patch-at -1 1))
  [
    move-to patch-at 0 1
    ask matrix-here [die]
  ]
  [ifelse (any? matrix-on patch-at 1 1 and (any? cells-on patch-at 0 1 or any? cells-on patch-at 1 0))
    [
      move-to patch-at 1 1
      set ycor ycor - 0.5
      ask matrix-here [die]
    ]
    [ifelse (any? matrix-on patch-at 1 0 and (any? cells-on patch-at 1 1 or any? cells-on patch-at 0 -1))
      [
        move-to patch-at 1 0
        set ycor ycor - 0.5
        ask matrix-here [die]
      ]
      [ifelse (any? matrix-on patch-at 0 -1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at -1 0))
        [
          move-to patch-at 0 -1
          ask matrix-here [die]
        ]
        [ifelse (any? matrix-on patch-at -1 0 and (any? cells-on patch-at 0 -1 or any? cells-on patch-at -1 1))
          [
            move-to patch-at -1 0
            set ycor ycor - 0.5
            ask matrix-here [die]
          ]
          [if (any? matrix-on patch-at -1 1 and (any? cells-on patch-at -1 0 or any? cells-on patch-at 0 1))
            [
              move-to patch-at -1 1
              set ycor ycor - 0.5
              ask matrix-here [die]
            ]
          ]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    ;; update neighborhood for cell
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule6-odd
  ifelse (not any? turtles-on patch-at 0 1 and (any? cells-on patch-at 1 1 or any? cells-on patch-at -1 1))
  [
    move-to patch-at 0 1
  ]
  [ifelse (not any? turtles-on patch-at 1 1 and (any? cells-on patch-at 0 1 or any? cells-on patch-at 1 0))
    [
      move-to patch-at 1 1
      set ycor ycor - 0.5
    ]
    [ifelse (not any? turtles-on patch-at 1 0 and (any? cells-on patch-at 1 1 or any? cells-on patch-at 0 -1))
      [
        move-to patch-at 1 0
        set ycor ycor - 0.5
      ]
      [ifelse (not any? turtles-on patch-at 0 -1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at -1 0))
        [
          move-to patch-at 0 -1
        ]
        [ifelse (not any? turtles-on patch-at -1 0 and (any? cells-on patch-at 0 -1 or any? cells-on patch-at -1 1))
          [
            move-to patch-at -1 0
            set ycor ycor - 0.5
          ]
          [if (not any? turtles-on patch-at -1 1 and (any? cells-on patch-at -1 0 or any? cells-on patch-at 0 1))
            [
              move-to patch-at -1 1
              set ycor ycor - 0.5
            ]
          ]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    ;; update neighborhood for cell
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule6-even ;; if there is only one cell neighbor and no matrix, move next to neighbor cell
  ifelse (not any? turtles-on patch-at 0 1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at -1 0))
  [
    move-to patch-at 0 1
    set ycor ycor - 0.5
  ]
  [ifelse (not any? turtles-on patch-at 1 0 and (any? cells-on patch-at 0 1 or any? cells-on patch-at 1 -1))
    [
      move-to patch-at 1 0
    ]
    [ifelse (not any? turtles-on patch-at 1 -1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at 0 -1))
      [
        move-to patch-at 1 -1
      ]
      [ifelse (not any? turtles-on patch-at 0 -1 and (any? cells-on patch-at 1 -1 or any? cells-on patch-at -1 -1))
        [
          move-to patch-at 0 -1
          set ycor ycor - 0.5
        ]
        [ifelse (not any? turtles-on patch-at -1 -1 and (any? cells-on patch-at 0 -1 or any? cells-on patch-at -1 0))
          [
            move-to patch-at -1 -1
          ]
          [if (not any? turtles-on patch-at -1 0 and (any? cells-on patch-at -1 -1 or any? cells-on patch-at 0 1))
            [
              move-to patch-at -1 0
            ]
          ]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    ;; update neighborhood for cell
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule7 ;; if cell is connected but only has matrix and empty neighbors, move to matrix with cell neighbors
  ifelse any? matrix-on neighbors-6 with [count cells in-radius 1.5 >= 2]
  [
    move-to one-of matrix-on neighbors-6 with [count cells in-radius 1.5 >= 2]
    ask matrix-here [die]
  ]
  [
    ifelse any? matrix-on neighbors-6 with [count cells in-radius 2 >= 2]
    [
      move-to one-of matrix-on neighbors-6 with [count cells in-radius 2 >= 2]
      ask matrix-here [die]
    ]
    [
      move-to one-of matrix-on neighbors-6
      ask matrix-here [die]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule8 ;; when a cell has two or less cell neighbors and only one matrix neighbor, it will move to the matrix neighbor
  move-to one-of matrix-on neighbors-6
  ask matrix-here [die]
  ifelse pxcor mod 2 = 0 [
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule9 ;; if no matrix neighbors and only one cel neighbor, move to empty patch in neighborhood with the most cells around
  ifelse any? patches in-radius 1.5 with [count cells in-radius 1.5 = 6 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
[
  move-to one-of patches in-radius 1.5 with [count cells in-radius 1.5 = 6 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
  if pxcor mod 2 = 0 [set ycor ycor - 0.5]
  ]
  [
    ifelse any? patches in-radius 1.5 with [count cells in-radius 1.5 = 5 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
    [
      move-to one-of patches in-radius 1.5 with [count cells in-radius 1.5 = 5 and (not any? turtles-here) and (not any? turtles-on patch-at 0 -0.5)]
      if pxcor mod 2 = 0 [set ycor ycor - 0.5]
    ]
    [
      ifelse any? patches in-radius 1.5 with [count cells in-radius 1.5 = 4 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
      [
        move-to one-of patches in-radius 1.5 with [count cells in-radius 1.5 = 4 and (not any? turtles-here) and (not any? turtles-on patch-at 0 -0.5)]
        if pxcor mod 2 = 0 [set ycor ycor - 0.5]
      ]
      [
        ifelse any? patches in-radius 1.5 with [count cells in-radius 1.5 = 3 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
        [
          move-to one-of patches in-radius 1.5 with [count cells in-radius 1.5 = 3 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
          if pxcor mod 2 = 0 [set ycor ycor - 0.5]
        ]
        [
          ifelse any? patches in-radius 1.5 with [count cells in-radius 1.5 = 2 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
          [
            move-to one-of patches in-radius 1.5 with [count cells in-radius 1.5 = 2 and not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
            if pxcor mod 2 = 0 [set ycor ycor - 0.5]
          ]
          [
            move-to one-of patches in-radius 1.5 with [not any? turtles-here and not any? turtles-on patch-at 0 -0.5]
            if pxcor mod 2 = 0 [set ycor ycor - 0.5]
          ]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
    set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
  ][
    set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to rule10-odd ;; cell neighbors = 1 , matrix = 4, and 1 empty space - move to empty space
  ifelse not any? turtles-on patch-at 0 1
    [move-to patch-at 0 1]
  [
    ifelse not any? turtles-on patch-at 1 1
    [move-to patch-at 1 1
      set ycor ycor - 0.5
    ]
    [
      ifelse not any? turtles-on patch-at 1 0
      [move-to patch-at 1 0
        set ycor ycor - 0.5
      ]
      [
        ifelse not any? turtles-on patch-at 0 -1
        [move-to patch-at 0 -1]
        [
          ifelse not any? turtles-on patch-at -1 0
          [move-to patch-at -1 0
            set ycor ycor - 0.5
          ]
          [
            if not any? turtles-on patch-at -1 1
            [move-to patch-at -1 1
              set ycor ycor - 0.5
            ]
          ]
        ]
      ]
    ]
  ]
end

to rule10-even ;; cell neighbors = 1 , matrix = 4, and 1 empty space - move to empty space
  ifelse not any? turtles-on patch-at 0 1
    [move-to patch-at 0 1
      set ycor ycor - 0.5
  ]
  [
    ifelse not any? turtles-on patch-at 1 0
    [move-to patch-at 1 0]
    [
      ifelse not any? turtles-on patch-at 1 -1
      [move-to patch-at 1 -1]
      [
        ifelse not any? turtles-on patch-at 0 -1
        [move-to patch-at 0 -1
          set ycor ycor - 0.5
        ]
        [
          ifelse not any? turtles-on patch-at -1 -1
          [move-to patch-at -1 -1]
          [
            if not any? turtles-on patch-at -1 0
            [move-to patch-at -1 0]
          ]
        ]
      ]
    ]
  ]
end

to assess-state-even
  if (pxcor mod 2 = 0) and (count cells-here = 1) and (pycor > min-pycor) and (pycor < max-pycor) and (pxcor > min-pxcor) and (pxcor < max-pxcor)
  [
    ifelse not any? cells-on neighbors-6
    [ifelse not any? matrix-on neighbors-6
      [die]
      [ifelse count matrix-on neighbors-6 = 6
        [hatch-cells 1
          [create-cell-link-with myself [set color blue]
           axiom3-divide-replacematrix
           set cycles-matrix-contact 0]
        ]
        [if (count matrix-on neighbors-6) <= 4 and (count cells-on neighbors-6 = 0)
          [
            axiom8-axiom7-divide]
          ]
      ]
    ]
    [ifelse count cells-on neighbors-6 = 6
      [die]
      [ifelse not any? matrix-on neighbors-6
        [ifelse count cells-on neighbors-6 = 1
          [
           hatch-matrix 1 [ set color white
            create-cell-link-with myself [set color blue]]
           axiom4-move-matrixbetween
           set cycles-matrix-contact 0
          ]
          ;; axiom 5:
          [
            die
          ]
        ]
        [ifelse (count matrix-on neighbors-6 + count cells-on neighbors-6) = 6 and (count matrix-on neighbors-6 < 6)
          [hatch-cells 1
           [create-cell-link-with myself [set color blue]
            axiom6-move
            set cycles-matrix-contact 0]
          ]
          [ifelse (count matrix-on neighbors-6 + count cells-on neighbors-6) <= 4
            [check-axiom8
             ;;axiom8-axiom7-divide
            ]
            [set still_count still_count + 1
              if still_count > 100 [set ss? True]]
          ]
        ]
      ]
    ]
    ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]
end

to assess-state-odd
  if (pxcor mod 2 != 0) and (count cells-here = 1) and (pycor > min-pycor) and (pycor < max-pycor) and (pxcor > min-pxcor) and (pxcor < max-pxcor)
  [
    ifelse not any? cells-on neighbors-6
    [ifelse not any? matrix-on neighbors-6
      [die]
      [ifelse count matrix-on neighbors-6 = 6
        [hatch-cells 1
          [create-cell-link-with myself [set color blue]
           axiom3-divide-replacematrix
          set cycles-matrix-contact 0]
        ]
        [if (count matrix-on neighbors-6) <= 4 and (count cells-on neighbors-6 = 0)
          [
            axiom8-axiom7-divide-odd]
          ]
      ]
    ]
    [ifelse count cells-on neighbors-6 = 6
      [die]
      [ifelse not any? matrix-on neighbors-6
        [ifelse count cells-on neighbors-6 = 1
          [
           hatch-matrix 1 [ set color white
            create-cell-link-with myself [set color blue]]
           axiom4-move-matrixbetween-odd
            set cycles-matrix-contact 0
          ]
          ;; axiom 5:
          [
            die
          ]
        ]
        [ifelse (count matrix-on neighbors-6 + count cells-on neighbors-6) = 6 and (count matrix-on neighbors-6 < 6)
          [hatch-cells 1
           [create-cell-link-with myself [set color blue]
            axiom6-move-odd
            set cycles-matrix-contact 0]
          ]
          [ifelse (count matrix-on neighbors-6 + count cells-on neighbors-6) <= 4
            [check-axiom8-odd
             ;;axiom8-axiom7-divide-odd
            ]
            [set still_count still_count + 1
              if still_count > 100 [set ss? True]]
          ]
        ]
      ]
    ]
    ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]
end

to axiom3-divide-replacematrix
    let matrix-patches neighbors-6 with [any? matrix-here]
    if any? matrix-patches
    [let target one-of matrix-patches
      face target
      move-to target
      ask target [die]
    ]
  ifelse pxcor mod 2 = 0 [
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]

end



to axiom4-move-matrixbetween
  ifelse any? cells-on patch-at 0 1
  [move-to patch-at 0 -1
    set ycor ycor - 0.5
  ]
  [ifelse any? cells-on patch-at 1 0
    [move-to patch-at -1 -1]
    [ifelse any? cells-on patch-at 1 -1
      [move-to patch-at -1 0]
      [ifelse any? cells-on patch-at 0 -1
        [move-to patch-at 0 1
          set ycor ycor - 0.5
        ]
        [ifelse any? cells-on patch-at -1 -1
          [move-to patch-at 1 0]
          [move-to patch-at 1 -1]
        ]
      ]
    ]
  ]
  ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to axiom4-move-matrixbetween-odd
  ifelse any? cells-on patch-at 0 1
  [move-to patch-at 0 -1
  ]
  [ifelse any? cells-on patch-at 1 1
    [move-to patch-at -1 0
    set ycor ycor - 0.5]
    [ifelse any? cells-on patch-at 1 0
      [move-to patch-at -1 1
      set ycor ycor - 0.5]
      [ifelse any? cells-on patch-at 0 -1
        [move-to patch-at 0 1
        ]
        [ifelse any? cells-on patch-at -1 0
          [move-to patch-at 1 1
          set ycor ycor - 0.5]
          [move-to patch-at 1 0
          set ycor ycor - 0.5]
        ]
      ]
    ]
  ]
 ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to axiom8-axiom7-divide
hatch-cells 1 [
create-cell-link-with myself [set color blue]
set cycles-matrix-contact 0
ifelse (not any? turtles-on patch-at 0 1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at -1 0))
[move-to patch-at 0 1
 set ycor ycor - 0.5
 ]
[ifelse (not any? turtles-on patch-at 1 0 and (any? matrix-on patch-at 0 1 or any? matrix-on patch-at 1 -1))
  [move-to patch-at 1 0
   ]
  [ifelse (not any? turtles-on patch-at 1 -1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at 0 -1))
    [move-to patch-at 1 -1
     ]
    [ifelse (not any? turtles-on patch-at 0 -1 and (any? matrix-on patch-at 1 -1 or any? matrix-on patch-at -1 -1))
      [move-to patch-at 0 -1
       set ycor ycor - 0.5
       ]
      [ifelse (not any? turtles-on patch-at -1 -1 and (any? matrix-on patch-at 0 -1 or any? matrix-on patch-at -1 0))
        [move-to patch-at -1 -1
         ]
        [if (not any? turtles-on patch-at -1 0 and (any? matrix-on patch-at -1 -1 or any? matrix-on patch-at 0 1))
          [move-to patch-at -1 0
           ]
          ]
      ]
    ]
  ]
]
 ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]

end

to axiom8-axiom7-divide-odd
hatch-cells 1 [
create-cell-link-with myself [set color blue]
set cycles-matrix-contact 0
ifelse (not any? turtles-on patch-at 0 1 and (any? matrix-on patch-at 1 1 or any? matrix-on patch-at -1 1))
[move-to patch-at 0 1
 ]
[ifelse (not any? turtles-on patch-at 1 1 and (any? matrix-on patch-at 0 1 or any? matrix-on patch-at 1 0))
  [move-to patch-at 1 1
       set ycor ycor - 0.5
   ]
  [ifelse (not any? turtles-on patch-at 1 0 and (any? matrix-on patch-at 1 1 or any? matrix-on patch-at 0 -1))
    [move-to patch-at 1 0
         set ycor ycor - 0.5
     ]
    [ifelse (not any? turtles-on patch-at 0 -1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at -1 0))
      [move-to patch-at 0 -1
       ]
      [ifelse (not any? turtles-on patch-at -1 0 and (any? matrix-on patch-at 0 -1 or any? matrix-on patch-at -1 1))
        [move-to patch-at -1 0
         set ycor ycor - 0.5
         ]
        [if (not any? turtles-on patch-at -1 1 and (any? matrix-on patch-at -1 0 or any? matrix-on patch-at 0 1))
          [move-to patch-at -1 1
           set ycor ycor - 0.5
           ]
          ]
      ]
    ]
  ]
]
    ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
]
end

to axiom6-move
ifelse (not any? cells-on patch-at 0 1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at -1 0))
[move-to patch-at 0 1
    set ycor ycor - 0.5
 ask matrix-here [die]]
[ifelse (not any? cells-on patch-at 1 0 and (any? cells-on patch-at 0 1 or any? cells-on patch-at 1 -1))
  [move-to patch-at 1 0
   ask matrix-here [die]]
  [ifelse (not any? cells-on patch-at 1 0 and (any? cells-on patch-at 1 0 or any? cells-on patch-at 0 -1))
    [move-to patch-at 1 0
     ask matrix-here [die]]
    [ifelse (not any? cells-on patch-at 0 -1 and (any? cells-on patch-at 1 -1 or any? cells-on patch-at -1 -1))
      [move-to patch-at 0 -1
          set ycor ycor - 0.5
       ask matrix-here [die]]
      [ifelse (not any? cells-on patch-at -1 -1 and (any? cells-on patch-at 0 -1 or any? cells-on patch-at -1 0))
        [move-to patch-at -1 -1
         ask matrix-here [die]]
        [if (not any? cells-on patch-at -1 0 and (any? cells-on patch-at -1 -1 or any? cells-on patch-at 0 1))
          [move-to patch-at -1 0
           ask matrix-here [die]]
          ]
      ]
    ]
  ]
]
   ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to axiom6-move-odd
ifelse (not any? cells-on patch-at 0 1 and (any? cells-on patch-at 1 1 or any? cells-on patch-at -1 1))
[move-to patch-at 0 1
 ask matrix-here [die]]
[ifelse (not any? cells-on patch-at 1 1 and (any? cells-on patch-at 0 1 or any? cells-on patch-at 1 0))
  [move-to patch-at 1 1
      set ycor ycor - 0.5
   ask matrix-here [die]]
  [ifelse (not any? cells-on patch-at 1 0 and (any? cells-on patch-at 1 1 or any? cells-on patch-at 0 -1))
    [move-to patch-at 1 0
        set ycor ycor - 0.5
     ask matrix-here [die]]
    [ifelse (not any? cells-on patch-at 0 -1 and (any? cells-on patch-at 1 0 or any? cells-on patch-at -1 0))
      [move-to patch-at 0 -1
       ask matrix-here [die]]
      [ifelse (not any? cells-on patch-at -1 0 and (any? cells-on patch-at 0 -1 or any? cells-on patch-at -1 1))
        [move-to patch-at -1 0
            set ycor ycor - 0.5
         ask matrix-here [die]]
        [if (not any? cells-on patch-at -1 1 and (any? cells-on patch-at -1 0 or any? cells-on patch-at 0 1))
          [move-to patch-at -1 1
           set ycor ycor - 0.5
           ask matrix-here [die]]
          ]
      ]
    ]
  ]
]
   ifelse pxcor mod 2 = 0 [
            ;; update neighborhood for cell
        set neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
end

to check-axiom8
ifelse (not any? turtles-on patch-at 0 1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at -1 0))
[axiom8-axiom7-divide]
[ifelse (not any? turtles-on patch-at 1 0 and (any? matrix-on patch-at 0 1 or any? matrix-on patch-at 1 -1))
  [axiom8-axiom7-divide]
  [ifelse (not any? turtles-on patch-at 1 -1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at 0 -1))
    [axiom8-axiom7-divide]
    [ifelse (not any? turtles-on patch-at 0 -1 and (any? matrix-on patch-at 1 -1 or any? matrix-on patch-at -1 -1))
      [axiom8-axiom7-divide]
      [ifelse (not any? turtles-on patch-at -1 -1 and (any? matrix-on patch-at 0 -1 or any? matrix-on patch-at -1 0))
        [axiom8-axiom7-divide]
        [ifelse (not any? turtles-on patch-at -1 0 and (any? matrix-on patch-at -1 -1 or any? matrix-on patch-at 0 1))
          [axiom8-axiom7-divide]
          [set still_count still_count + 1
              if still_count > 20 [set ss? True]]
          ]
      ]
    ]
  ]
  ]
end

to check-axiom8-odd
ifelse (not any? turtles-on patch-at 0 1 and (any? matrix-on patch-at 1 1 or any? matrix-on patch-at -1 1))
[axiom8-axiom7-divide-odd]
[ifelse (not any? turtles-on patch-at 1 1 and (any? matrix-on patch-at 0 1 or any? matrix-on patch-at 1 0))
  [axiom8-axiom7-divide-odd]
  [ifelse (not any? turtles-on patch-at 1 0 and (any? matrix-on patch-at 1 1 or any? matrix-on patch-at 0 -1))
    [axiom8-axiom7-divide-odd]
    [ifelse (not any? turtles-on patch-at 0 -1 and (any? matrix-on patch-at 1 0 or any? matrix-on patch-at -1 0))
      [axiom8-axiom7-divide-odd]
      [ifelse (not any? turtles-on patch-at -1 0 and (any? matrix-on patch-at 0 -1 or any? matrix-on patch-at -1 1))
        [axiom8-axiom7-divide-odd]
        [ifelse (not any? turtles-on patch-at -1 1 and (any? matrix-on patch-at -1 0 or any? matrix-on patch-at 0 1))
          [axiom8-axiom7-divide-odd]
            [set still_count still_count + 1
              if still_count > 20 [set ss? True]]
          ]
      ]
    ]
  ]
]
end


to find-cysts
  ask patches with [not any? turtles-here]
  [ifelse count cells-on patch-neighbors-6 >= 1
    [set pcolor blue]
    [let empty-patches patch-neighbors-6 with [pcolor = blue]
      if empty-patches != nobody
      [set pcolor blue]
  ]
  ]
end

to find-lumen
  ask patches with [pcolor = blue and count lumen-here = 0]
  [
   sprout-lumen 1
    [set color yellow
     set shape "circle"
     set group-id 0
     set unvisited? True
     if pxcor mod 2 = 0  [ set ycor ycor - 0.5 ]
     ;; define neighborhood of patches
     ifelse pxcor mod 2 = 0 [
        set lumen-neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set lumen-neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]
  ]

  ask lumen [
    ifelse pxcor mod 2 = 0 [
        set lumen-neighbors-6 turtles-on patches at-points [[0 1] [1 0] [1 -1] [0 -1] [-1 -1] [-1 0]]
    ][
        set lumen-neighbors-6 turtles-on patches at-points [[0 1] [1 1] [1  0] [0 -1] [-1  0] [-1 1]] ]
  ]
  ask lumen [
    create-lumen-links-with lumen-on lumen-neighbors-6 with [not link-neighbor? myself]
    [set color orange
      set thickness 0.3
    ]
  ]

;  ask lumen [
;    if count cells-on lumen-neighbors-6 >= 1 [
;    let cyst-cells cells-on lumen-neighbors-6
;      if cyst-cells != nobody
;      [ask cyst-cells [
;        if counted? = False [set counted? True]
;        ]
;      ]
;    ]
;  ]
end

to walk-up
  while [count lumen-on patch-at 0 1 = 1 ] [
    move-to patch-at 0 1
    if pxcor mod 2 = 0  [ set ycor ycor - 0.5 ]
    set up-moves up-moves + 1
    ask lumen-here [set unvisited? False]
    if any? lumen-here with [group-id > 0]
    [set id-in-proximity? True]
    if any? lumen in-radius 1.5 with [group-id > 0]
    [
      set id-in-proximity? True
    ]
    ask lumen-here
    [if (count lumen in-radius 1.5 with [xcor != [xcor] of myself and unvisited? = True] >= 1)
      [
        set gateway? True
        set gateway-visited? False
      ]
    ]
  ]
end

to walk-down
  while [count lumen-on patch-at 0 -1 = 1 ] [
    move-to patch-at 0 -1
    if pxcor mod 2 = 0  [ set ycor ycor - 0.5 ]
    set down-moves down-moves + 1
    if any? lumen-here with [group-id > 0]
    [set id-in-proximity? True]
  ]
end

to walk-back-down
  while [up-moves > 0]
  [
    move-to patch-at 0 -1
    if pxcor mod 2 = 0  [ set ycor ycor - 0.5 ]
    set up-moves up-moves - 1
  ]
end

to walk-back-up
  while [down-moves > 0]
  [
    move-to patch-at 0 1
    if pxcor mod 2 = 0  [ set ycor ycor - 0.5 ]
    set down-moves down-moves - 1
  ]
end


to walk-up-and-id
  while [count lumen-on patch-at 0 1 = 1 ] [
    move-to patch-at 0 1
    if pxcor mod 2 = 0  [ set ycor ycor - 0.5 ]
    if any? lumen-here ;with [group-id = 0]
    [ask lumen-here
      [
        set group-id [group-id] of myself
        set color pink
        ask lumen-on lumen-neighbors-6
        [
          set group-id [group-id] of myself
          set color pink
        ]
      ]
    ]
  ]
end

to walk-down-and-id
  while [count lumen-on patch-at 0 -1 = 1 ] [
    move-to patch-at 0 -1
    if pxcor mod 2 = 0  [ set ycor ycor - 0.5 ]
    if any? lumen-here ;with [group-id = 0]
    [ask lumen-here
      [
        set group-id [group-id] of myself
        set color pink
        ask lumen-on lumen-neighbors-6
        [
          set group-id [group-id] of myself
          set color pink
        ]
      ]
    ]
  ]
end

to make-cyst-id
set grid-x-pos min-pxcor
  while [grid-x-pos <= max-pxcor] [
    set grid-y-pos max-pycor
    while [grid-y-pos >= min-pycor] [
      ask patch grid-x-pos grid-y-pos
      [
        if any? lumen-here
        [
          ask lumen-here
          [
            ifelse group-id = 0
            [
              hatch-walker 1
              [
                set id-in-proximity? False
                walk-up
                walk-back-down
                walk-down
                walk-back-up
              ]
              if any? walker-here with [id-in-proximity? = False]
              [
                set group-id num-cysts
                set num-cysts num-cysts + 1
                set color pink
                ask lumen-on lumen-neighbors-6
                [
                  set group-id [group-id] of myself
                  set color pink
                ]
                ask walker-here
                [
                  set group-id [group-id] of myself
                  walk-up-and-id
                  walk-down-and-id
                  die
                ]
              ]
            ]
            [
              if group-id > 0
              [
                ask lumen-on lumen-neighbors-6
                [
                  set group-id [group-id] of myself
                  set color pink
                ]
                hatch-walker 1
                [
                  set group-id [group-id] of myself
                  walk-up-and-id
                  walk-down-and-id
                  die
                ]
              ]
            ]
          ]
        ]
      ]
      set grid-y-pos grid-y-pos - 1
    ]
    set grid-x-pos grid-x-pos + 1
  ]

end

to make-cell-id
  ask lumen
  [
    ask cells-on lumen-neighbors-6
    [
      set group-id [group-id] of myself
    ]
  ]
end

to propagate-id
  ask cells
  [
    ask cells-on neighbors-6
    [
      set group-id [group-id] of myself
    ]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
5
30
837
880
-1
-1
8.41
1
10
1
1
1
0
0
0
1
-48
49
-49
50
1
1
1
ticks
60.0

BUTTON
885
495
1007
535
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

TEXTBOX
885
476
1024
494
Run model
11
0.0
0

CHOOSER
885
405
1010
450
atom-shape
atom-shape
"hex" "hexline" "thin-line" "line" "spikes90" "default"
0

BUTTON
1015
405
1154
450
apply shape
ask cells [set shape atom-shape]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1012
495
1147
535
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
890
185
970
268
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
885
370
1175
388
_____________________________________________
11
0.0
1

TEXTBOX
885
455
1180
481
_____________________________________________
11
0.0
1

TEXTBOX
857
105
882
141
1
30
14.0
1

TEXTBOX
858
385
878
421
2
30
14.0
1

TEXTBOX
858
485
883
521
6
30
14.0
1

SLIDER
885
145
1082
178
overlay-matrix-percentage
overlay-matrix-percentage
0
100
100.0
1
1
NIL
HORIZONTAL

SLIDER
1180
205
1352
238
num-matrix-diff
num-matrix-diff
1
6
4.0
1
1
NIL
HORIZONTAL

SWITCH
1180
165
1297
198
diff-matrix?
diff-matrix?
0
1
-1000

SLIDER
1180
250
1352
283
cycles-diff-matrix
cycles-diff-matrix
0
20
3.0
1
1
NIL
HORIZONTAL

SWITCH
1180
350
1312
383
diff-induction?
diff-induction?
0
1
-1000

SLIDER
1180
395
1352
428
num-diff-ind
num-diff-ind
0
6
2.0
1
1
NIL
HORIZONTAL

SLIDER
1180
440
1352
473
cycles-diff-ind
cycles-diff-ind
1
20
4.0
1
1
NIL
HORIZONTAL

TEXTBOX
1177
80
1282
116
3
30
14.0
1

TEXTBOX
1185
315
1335
351
4
30
14.0
1

TEXTBOX
1180
295
1420
321
__________________________________________
9
0.0
1

SLIDER
887
105
1059
138
num-cells
num-cells
0
400
97.0
1
1
NIL
HORIZONTAL

SLIDER
1177
540
1364
573
undiff-num-inhibition
undiff-num-inhibition
1
6
3.0
1
1
NIL
HORIZONTAL

BUTTON
887
575
979
608
GO-ONCE
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1182
505
1332
541
5
30
14.0
1

TEXTBOX
1182
490
1437
516
____________________________________
11
0.0
1

SLIDER
1182
120
1389
153
max-ticks-differentiation
max-ticks-differentiation
0
20
2.0
1
1
NIL
HORIZONTAL

CHOOSER
890
50
1028
95
Rule-set
Rule-set
"Original" "New"
1

SLIDER
890
280
1062
313
max-divisions
max-divisions
1
5
2.0
1
1
NIL
HORIZONTAL

CHOOSER
890
325
1028
370
culture-condition
culture-condition
"embedded" "clustered"
1

SLIDER
1115
35
1287
68
cluster-size
cluster-size
0.5
2
2.0
0.5
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

Most materials are not continuous arrangements of atoms. Rather, they are composed of thousands or millions of microscopic crystals, known as grains.  This model shows how the configuration and sizes of these grains change over time.  Grain size is a very important characteristic for evaluating the mechanical properties of materials; it is exhaustively studied in metallurgy and materials science.

Usually this kind of study is made by careful analysis and comparison of pictures taken in microscopes, sometimes with the help of image analysis software.  Recently, as the processing power of computers has increased, a new and promising approach has been made possible: computer simulation of grain growth.

Anderson, Srolovitz et al. proposed the most widely known and employed theory for computer modeling and simulation of grain growth, using the Monte Carlo method. Instead of considering the grains as spheres, and being obliged to make numerous geometrical approximations, Anderson proposed that the computer would simulate the behavior of each individual atom in the system. Each atom would follow a very simple rule: it will always try to have, in its immediate neighborhood, as many atoms as possible with the same orientation as it. It will do so by randomly (hence Monte Carlo) re-orienting itself and seeing if it is more stable than it was before. If it is, it will stay in its new orientation, and if not, it will revert back to its previous orientation.

This model is part of the MaterialSim (Blikstein & Wilensky, 2004) curricular package. To learn more about MaterialSim, see http://ccl.northwestern.edu/rp/materialsim/index.shtml.

## HOW IT WORKS

The basic algorithm of the simulation is simple: each atom continuously tries to be as stable as possible.  Its stability is based on the number of neighbors with similar orientations: the more similar neighbors, the more stable it is. If it has only few similar neighbors, it will try to relocate to a more stable position.  The steps for each tick in the model are:

1. Choose a random atom.
2. Ask that atom to calculate its present energy (based on its stability).  The atom does this by counting how many similar neighbors it has.
3. The atom then chooses at random one of its neighbors, and orients itself in the same direction as that neighbor.
4. The same atom then calculates its energy, based on this new, tentative orientation.
5. Finally the atom compares the energy levels in each of the two states: the lowest value "wins", i.e., the more similar neighbors, the more stable the atom is.
6. Repeat steps 1-6.

The **annealing-temperature** slider controls the probability of maintaining a reorientation that yields less stability.  The **percent-element2** slider defines the percentage of second-phase particles to be created when the user setups the simulation.  Those particles are not movable and are not subject to grain growth. Atoms with element2-particles as neighbors will see them as dissimilar.

Note that the actual number of atoms is small compared to a real metal sample.  Also, real materials are three-dimensional, while this model is 2D.

## HOW TO USE IT

### (1) Simulation starting point.
**starting-point**: You can start from a random arrangement or a picture from a microscope.  File formats accepted are: `.jpg`, `.png`, `.bmp`, and `.gif`. The image will be automatically resized to fit into the world, but maintaining its original aspect ratio. Note that the image MUST HAVE THE SAME ASPECT RATIO AS THE WORLD. In other words, if the world is square, the image should be square as well. Prior to importing the image, it is recommended to clean it up using an image editing software (increase contrast, remove noise).  Try to experiment various combinations of values for the view's size and the patch size to get the best results

**percent-element2**: You can also determine if you want a certain percentage of a second element, sometimes called grain refiner.

**setup**: Sets up the model as indicated.

### (2) Change the shape of atoms.

Different shapes might help visualize the atomic planes, or the proportion of different types of grains.

**atom-shape**: Choose which shape you want.

**apply-shape**: Re-draws all the atoms with the shape you chose.

### (3) Draw/edit grains (optional)

You can use this to "draw" new grain structures, or edit the grains at any point during the simulation.

**draw**: When this button is pressed down, you can change the color (orientation) of atoms in the model.

**brush-size**: This determines how many atoms’ orientation you change at a time.

**erase-all**:  This changes the orientation of all atoms.

**draw-color**: Here you can change the color (and orientation) that that you change atoms to, when you draw on them.

### (4) Run Model.

You can choose to have the model keep repeating the steps (as described above) or to run just one tick. Running just one tick will allow you to see more clearly what happens over time.

**annealing-temperature**: changes the probability of non-favorable orientation flips to happen.  A 10% value, for instance, means that 10% of non-favorable flips will be maintained.  This mimics the effect of higher temperatures.

### Grain size plot and calculations

**ticks-per-measurement **: to increase the model's speed, the user can choose not to calculate grain size at every time step.  If grain size is calculated at every ten time units (20, 30, 40 etc.), the performance is slightly increased.  This only affects the plot and the monitors, but not the actual simulation.

**measure grains now**: if the user feels there is too much time between each grain measurement, and wants to calculate it immediately, this button can be used.

### Plots and monitors

**Grain Size (log-log)**: Grain size vs. time, in a log-log scale.  Under normal conditions (**annealing-temperature** = 0 and **percent-element2** = 0), this plot should be a straight line with an angular coefficient of approximately 0.5.

**Grain Size**: grain size

**Log Time**: log of ticks in the simulation so far

**Log Grain Size**: Log of the grain size

**growth-exponent**: the angular coefficient of the **Grain Size (log-log)** plot.  This number should approach 0.5 with **annealing-temperature** = 0 and **percent-element2** = 0.

## THINGS TO NOTICE

When you setup with a random orientation and run the simulation, notice that the speed of growth decreases with time.  Toward the end of the simulation, you might see just two or three grains that fight with each other for along time.  One will eventually prevail, but this logarithmic decrease of speed is an important characteristic of grain growth.  That is why the **Grain Size (log-log) plot is a straight line in a "log-log" scale.
Notice also that if you draw two grains, one concave and one convex, their boundary will tend to be a straight line, if you let the simulation run long enough.  Every curved boundary is unstable because many atoms at its interface will have more different than equal neighbors.

## THINGS TO TRY

Increase the value of the **annealing-temperature** slider. What happens to the **Grain Size (log-log)** plot, and to the boundaries' shapes?

Try to increase the **percent-element2** slider to 5%.  Then choose “random arrangement” and click **setup**, and **go**. What happens to grain growth? Now try several values (1, 3, 5, 7, 9%), for instance.  What happens with the final grain size? What about the **Grain Size (log-log) plot and the **Growth Exponent**?

One advanced use of this model would be to get a digital picture of a real metallic sample, reduce noise and increase contrast with image editing programs, and load into this model using the **Import Image** feature in the starting-point chooser.  Don't forget to update the width and height of the view's size to accommodate the picture, and also to change the patch size in order to be able to see the whole sample.

## EXTENDING THE MODEL

This model assumes that the misorientation between two grains has no effect on their growth rates.  Two grains with a very similar crystallographic orientation have the same growth rate as grains whose orientations differ by a lot.  Try to take the angular misorientation into consideration.

When we insert second-phase particles, all of them have the same size.  Try to create a slider that changes the size of these particles.

## NETLOGO FEATURES

Rather than containing all of the code that updates variables used for plotting, the **Grain Size (log-log)** plot calls a procedure that does this. The reason is that there is quite a lot of code, and it would be difficult to work with inside the plot. Notice also that this code updates the global variables shown in the monitors to the right of the plot.

The model uses a hexagonal grid as opposed to the usual, square one. It also uses different shapes for different visualization purposes. Finally, it uses the `import-pcolors` primitive for image import.

## RELATED MODELS

Crystallization Basic
Crystallization Directed
Crystallization Moving

## CREDITS AND REFERENCES

This model is part of the MaterialSim (Blikstein & Wilensky, 2004) curricular package. To learn more about MaterialSim, see http://ccl.northwestern.edu/rp/materialsim/index.shtml.

Two papers describing the use of this model in education are:
Blikstein, P. & Wilensky, U. (2005) Less is More: Agent-Based Simulation as a Powerful Learning Tool in Materials Science.  The IV International Conference on Autonomous Agents and Multiagent Systems. Utrecht, Netherlands.

Blikstein, P. & Wilensky, U. (2004) MaterialSim: An agent-based simulation toolkit for Materials Science learning. (PDF, 1.5 MB) Proceedings of the International Conference on Engineering Education. Gainesville, Florida.

The core algorithm of the model was developed at the University of Sao Paulo and published in:  Blikstein, P. and Tschiptschin, A. P. Monte Carlo simulation of grain growth (II). Materials Research, Sao Carlos, 2 (3), p. 133-138, jul. 1999.

Available for download at: http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1516-14391999000300004&lng=en&nrm=iso&tlng=en. See also http://www.pmt.usp.br/paulob/montecarlo/ for more information (in Portuguese).

## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Blikstein, P. and Wilensky, U. (2005).  NetLogo MaterialSim Grain Growth model.  http://ccl.northwestern.edu/netlogo/models/MaterialSimGrainGrowth.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 2005 Uri Wilensky.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.

<!-- 2005 Cite: Blikstein, P. -->
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circ
true
0
Circle -7500403 true true 10 11 278

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

hex
false
0
Polygon -7500403 true true 0 150 75 30 225 30 300 150 225 270 75 270

hexline
true
0
Polygon -7500403 true true 0 150 75 30 225 30 300 150 225 270 75 270
Rectangle -1 true false 121 47 182 252

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Rectangle -7500403 true true 135 0 165 315

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

rectangle
true
0
Polygon -7500403 true true 67 36 67 262 235 262 235 35

spikes90
true
0
Circle -7500403 true true 61 62 177
Line -7500403 true 135 66 131 154
Rectangle -7500403 true true 135 -4 166 68
Rectangle -7500403 true true 142 132 219 134
Rectangle -7500403 true true 196 136 304 165
Rectangle -7500403 true true -9 135 68 166
Rectangle -7500403 true true 131 176 132 226
Rectangle -7500403 true true 135 225 165 313

square
false
0
Rectangle -7500403 true true 15 15 286 285
Line -1 false 6 6 293 6
Line -1 false 293 6 293 293
Line -1 false 293 292 8 292
Line -1 false 8 292 8 6

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

t
true
0
Rectangle -7500403 true true 46 47 256 75
Rectangle -7500403 true true 135 76 167 297

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

thin-line
true
0
Line -7500403 true 150 0 150 300

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
