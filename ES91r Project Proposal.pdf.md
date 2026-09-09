# ES91r Project Proposal.pdf

_extracted by pdftext_

--- page 1 of 3 ---
Student  Name:  Andrew  Rodriguez  PI  Name:  Dr.  Nicholas  Todd  Primary  Mentor  Name:  Dr.  Nicholas  Todd  (day  to  day  contact  Dr.  Bernie  Owusu-Yaw)   Title  of  Project:  Machine  Learning  Based  Acoustic  Feedback  Controller  for  FUS-BBB  Mediated  
AAV
 
Delivery
 
  Project  Motivation:  The  blood-brain  barrier  (BBB)  is  a  major  obstacle  for  the  delivery  of  novel  
gene
 
therapies,
 
such
 
as
 
Adeno-Associated
 
Viruses
 
(AAVs),
 
to
 
the
 
central
 
nervous
 
system.
 
Focused
 
Ultrasound
 
Blood-Brain
 
Barrier
 
(FUS-BBB)
 
opening
 
is
 
a
 
promising,
 
non-invasive
 
technique
 
that
 
uses
 
targeted
 
acoustic
 
waves
 
and
 
circulating
 
microbubbles
 
to
 
temporarily
 
increase
 
BBB
 
permeability.
 
However,
 
ensuring
 
both
 
the
 
safety
 
and
 
efficacy
 
of
 
this
 
procedure
 
requires
 
precise
 
control
 
over
 
the
 
ultrasound
 
dose.
 
Currently,
 
Dr.
 
Todd’s
 
lab
 
has
 
a
 
need
 
for
 
real
 
time
 
monitoring
 
to
 
prevent
 
under-dosing
 
(ineffective
 
delivery)
 
or
 
over-dosing
 
(tissue
 
damage).
 
As
 
the
 
microbubbles
 
cavitate
 
during
 
FUS,
 
they
 
generate
 
distinct
 
acoustic
 
emissions
 
which
 
can
 
be
 
measured
 
by
 
the
 
piezoelectric
 
crystals
 
of
 
an
 
ultrasound
 
system.
 
The
 
overarching
 
goal
 
of
 
this
 
research
 
is
 
to
 
design
 
a
 
dynamic
 
feedback
 
controller
 
capable
 
of
 
processing
 
these
 
acoustic
 
emissions
 
in
 
real
 
time
 
to
 
predict
 
and
 
control
 
the
 
quantity
 
of
 
AAV
 
delivered
 
to
 
the
 
brain.
  Detailed  Description  of  Your  Specific  Project:  During  this  13  week  independent  study,  my  project  will  focus  on  the  first  two  phases  of  
designing
 
this
 
acoustic
 
feedback
 
controller,
 
heavily
 
utilizing
 
principles
 
of
 
signal
 
processing,
 
machine
 
learning,
 
and
 
quantitative
 
biology.
  Specifically,  I  will  process  raw,  time  series  acoustic  emission  data  generated  from  previous  FUS  
experiments
 
(some
 
of
 
this
 
data
 
has
 
already
 
been
 
collected,
 
some
 
more
 
may
 
be
 
collected
 
this
 
summer
 
by
 
Dr.
 
Todd’s
 
group).
 
Using
 
Python,
 
I
 
will
 
design
 
feature
 
extraction
 
pipelines,
 
such
 
as
 
the
 
Fourier
 
transform
 
(which
 
shows
 
preliminary
 
promise
 
on
 
the
 
second
 
harmonic),
 
to
 
isolate
 
the
 
specific
 
acoustic
 
signatures
 
correlated
 
with
 
BBB
 
opening.
 
Additionally,
 
I
 
will
 
engage
 
in
 
some
 
hands-on
 
measurement
 
and
 
experimental
 
validation
 
by
 
assisting
 
with
 
in
 
vivo
 
FUS
 
experiments
 
in
 
mouse
 
models.
 
Following
 
the
 
physical
 
experiments,
 
I
 
will
 
perform
 
tissue
 
sectioning,
 
staining,
 
and
 
fluorescence
 
microscopy
 
to
 
quantitatively
 
measure
 
the
 
ground
 
truth
 
concentration
 
of
 
AAVs
 
successfully
 
delivered
 
to
 
the
 
target
 
tissue.
  By  the  end  of  the  term,  I  will  integrate  these  two  workflows  by  applying  machine  learning  models  
to
 
map
 
the
 
extracted
 
acoustic
 
features
 
(the
 
input)
 
to
 
the
 
physical
 
AAV
 
delivery
 
outcomes
 
(the
 
output).
 
Through
 
this
 
project,
 
I
 
will
 
gain
 
practical
 
experience
 
in
 
bridging
 
computational
 
modeling
 
with
 
wet-lab
 
biological
 
measurement,
 
applying
 
my
 
coursework
 
in
 
Computer
 
Science
 
and
 
Biomedical
 
Engineering
 
to
 
a
 
translational
 
medicine
 
problem.
  Final  Deliverable:  A  detailed  written  report,  formatted  as  a  draft  journal  manuscript,  detailing  the  feature  extraction  
methodology,
 
the
 
experimental
 
protocol,
 
and
 
the
 
preliminary
 
predictive
 
accuracy
 
of
 
the
 
baseline

--- page 2 of 3 ---
machine  learning  model.  Additionally,  dependent  upon  need  in  the  fall  as  well  as  model  
accuracy,
 
I
 
will
 
make
 
the
 
model
 
easily
 
accessible
 
and
 
deployable
 
to
 
Dr.
 
Todd’s
 
group,
  Preliminary  Weekly  Schedule   Tuesday:  1:00  PM  -  5:00  PM   Wednesday:  10:00  AM  -  12:00  PM   Thursday:  10:00  AM  -  2:00  PM    Week  1:  Complete  lab  safety  training  (BWH/HMS),  review  relevant  literature  on  FUS-BBB  and  
acoustic
 
cavitation,
 
and
 
finalize
 
computational
 
environment
 
setup.
  Week  2:  Familiarize  with  existing  raw  acoustic  emission  datasets;  begin  writing  initial  Python  
scripts
 
for
 
data
 
import
 
and
 
basic
 
filtering.
  Week  3 :  Develop  and  test  feature  extraction  pipelines  (e.g.,  Fourier  transforms)  on  existing  
acoustic
 
data.
  Week  4:  Refine  feature  extraction  algorithms  to  isolate  harmonic,  subharmonic,  and  broadband  
emission
 
signatures.
  Week  5:  Begin  assisting  with  in  vivo  FUS  experiments;  observe  and  learn  the  focused  
ultrasound
 
protocols.
  Week  6:  Continue  assisting  with  in  vivo  experiments;  begin  training  on  wet  lab  techniques  for  
post
 
experiment
 
tissue
 
handling.
  Week  7:  Perform  brain  sectioning  and  biological  staining  on  tissue  from  recent  FUS  
experiments
 
to
 
prepare
 
for
 
AAV
 
quantification.
  Week  8:  Utilize  fluorescence  microscopy  and  imaging  software  to  quantitatively  measure  
physical
 
AAV
 
delivery
 
in
 
stained
 
tissue
 
samples.
  Week  9:  Align  ground  truth  AAV  measurement  data  (from  Week  8)  with  the  corresponding  
acoustic
 
emission
 
data
 
logs.
  Week  10:  Begin  building  and  training  the  baseline  machine  learning  model  using  the  newly  
structured
 
dataset
 
(mapping
 
acoustic
 
features
 
to
 
AAV
 
delivery).
  Week  11:  Test  and  evaluate  model  accuracy;  perform  iterative  adjustments  to  the  selected  
features
 
or
 
model
 
parameters.

--- page 3 of 3 ---
Week  12:  Finalize  data  analysis;  synthesize  findings  into  visualizations  (graphs,  performance  
metrics)
 
for
 
the
 
final
 
report.
  Week  13 :  Draft,  review,  and  finalize  the  summative  written  report  for  submission  by  the  end  of  
Reading
 
Period.
