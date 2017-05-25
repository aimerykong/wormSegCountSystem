This is the system for automated C. elegans analysis, including semantic segmentation, instance segmentation and counting.


Directory 'step1_annotation' is the interface to annotated worms. It contains two ways, manual annotation and semi-automatical annotation. In manual version, the users are supposed to draw a line along the worm body to annotate; while in the semi-automatical version, the users should click the worm head and tail (order does not matter), then a line is automatically drawn (chain model by dynamic programming). It may fail, so the users will be asked to double check before submitting the annotation.


Shu Kong @ UCI
05/24/2017
