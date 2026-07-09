## 11c. Validated composite score definitions

Each composite is a **prorated sum**: mean(available item scores) × n_items, requiring valid answers for more than half of items. Reverse-worded items (flagged per scale) are flipped on their own min/max before summing. Items are matched to survey responses by question text and merged across survey administrations. PHQ-9 and GAD-7 pool EHHWB and COPE administrations with EHHWB priority and a `from_cope` covariate; PSS-10 pools SDOH and COPE administrations with SDOH priority and the same source covariate. The score is then inverse-normal-transformed and residualized like any quantitative trait (§4.1). Phenotype ids are prefixed `comp_`.

### GAD-7 — Generalized Anxiety Disorder scale (anxiety)

- **Items:** 7
- **Per-item scoring:** Not at all = 0, Several days = 1, Over half the days = 2, Nearly all days = 3
- **Total score:** prorated sum of 7 items; no reverse-keyed items
- **Auto-built:** yes (comp_gad7_anxiety)
- **Questions:**
    - Feeling nervous, anxious, or on edge
    - Not being able to stop or control worrying
    - Worrying too much about different things
    - Trouble relaxing
    - Being so restless that it's hard to sit still
    - Becoming easily annoyed or irritable
    - Feeling afraid as if something awful might happen

### PHQ-9 — Patient Health Questionnaire (depression)

- **Items:** 9
- **Per-item scoring:** Not at all = 0, Several days = 1, Over half the days = 2, Nearly all days = 3
- **Total score:** prorated sum of 9 items; no reverse-keyed items
- **Auto-built:** yes (comp_phq9_depression)
- **Questions:**
    - Little interest or pleasure in doing things
    - Feeling down, depressed, or hopeless
    - Trouble falling or staying asleep, or sleeping too much
    - Feeling tired or having little energy
    - Poor appetite or overeating
    - Feeling bad about yourself or that you are a failure or have let yourself or your family down
    - Trouble concentrating on things, such as reading the newspaper or watching television
    - Moving or speaking so slowly that other people could have noticed? Or the opposite - being so fidgety or restless that you have been moving around a lot more than usual
    - Thoughts that you would be better off dead or of hurting yourself in some way

### PSS — Perceived Stress Scale

- **Items:** 10
- **Per-item scoring:** 2 answer scales across items (shown per item below)
- **Total score:** prorated sum of 10 items; 4 reverse-keyed
- **Auto-built:** yes (comp_pss_perceived_stress)
- **Pooling:** SDOH is primary; COPE fills COPE-only responses. The GWAS residualization includes `from_cope`.
- **Questions:**
    - In the last month, how often have you been upset because of something that happened unexpectedly?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt that you were unable to control the important things in your life?  — [Never=0.0, Almost never=1.0, Sometime=2.0, Fairly often=3.0, Often=4.0]
    - In the last month, how often have you felt nervous and "stressed?"  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt confident about your ability to handle your personal problems? *(reverse-keyed)*  — [Never=0.0, Almost never=1.0, Sometime=2.0, Fairly often=3.0, Often=4.0]
    - In the last month, how often have you felt that things were going your way? *(reverse-keyed)*  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you found that you could not cope with all the things that you had to do?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you been able to control irritations in your life? *(reverse-keyed)*  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt that you were on top of things? *(reverse-keyed)*  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you been angered because of things that were outside of your control?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]
    - In the last month, how often have you felt difficulties were piling up so high that you could not overcome them?  — [Never=0.0, Almost Never=1.0, Sometime=2.0, Fairly Often=3.0, Often=4.0, Sometimes=2.0, Very Often=4.0]

### ACE — Adverse Childhood Experiences

- **Items:** 11
- **Per-item scoring:** 3 answer scales across items (shown per item below)
- **Total score:** prorated sum of 11 items; no reverse-keyed items
- **Auto-built:** yes (comp_ace_adversity)
- **Questions:**
    - During your first 18 years of life, did you live with anyone who was depressed, mentally ill, or suicidal? (ACE category: Mentally ill household member)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, did you live with anyone who was a problem drinker or alcoholic? (ACE category: Substance abuse in household)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, did you live with anyone who used illegal street drugs or who abused prescription medications? (ACE category: Substance abuse in household)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, did you live with anyone who served time or was sentenced to serve time in a prison, jail, or other correctional facility? (ACE category: Incarcerated household member)  — [Yes=1.0, No=0.0]
    - During your first 18 years of life, were your parents separated or divorced? (ACE category: Parental separation/divorce)  — [Yes=1.0, No=0.0, Parents not married=0.0]
    - During your first 18 years of life, how often did your parents or adults in your home ever slap, hit, kick, punch or beat each other up? (ACE category: Violence between adults in household)  — [Never=0.0, Once=1.0, More than once=1.0]
    - Before age 18, how often did a parent or adult in your home ever hit, beat, kick, or physically hurt you in any way? Do not include spanking. (ACE category: Physical abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did a parent or adult in your home ever swear at you, insult you, or put you down? (ACE category: Emotional abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did anyone at least 5 years older than you or an adult, ever touch you sexually? (ACE category: Sexual abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did anyone at least 5 years older than you or an adult, try to make you touch them sexually? (ACE category: Sexual abuse)  — [Never=0.0, Once=1.0, More than once=1.0]
    - During your first 18 years of life, how often did anyone at least 5 years older than you or an adult, force you to have sex? (ACE category: Sexual abuse)  — [Never=0.0, Once=1.0, More than once=1.0]

### IES — Impact of Event Scale (event-related distress)

- **Items:** 6
- **Per-item scoring:** Not at all = 0, A little bit = 1, Moderately = 2, Quite a bit = 3, Extremely = 4
- **Total score:** prorated sum of 6 items; no reverse-keyed items
- **Auto-built:** yes (comp_ies_event_impact)
- **Questions:**
    - In the past 7 days, I thought about COVID-19 when I didn't mean to.
    - In the past 7 days, I felt watchful or on-guard.
    - In the past 7 days, other things kept making me think about COVID-19.
    - In the past 7 days, I was aware that I still had a lot of feelings about COVID-19, but I didn't deal with them.
    - In the past 7 days, I tried not to think about COVID-19.
    - In the past 7 days, I had trouble concentrating.

### ASRS — Adult ADHD Self-Report Scale (Part A screener)

- **Items:** 6
- **Per-item scoring:** 2 answer scales across items (shown per item below)
- **Total score:** prorated sum of 6 items; no reverse-keyed items
- **Auto-built:** yes (comp_asrs_adhd)
- **Questions:**
    - How often do you have trouble wrapping up the final details of a project, once the challenging parts have been done?  — [Never=0.0, Rarely=0.0, Sometimes=1.0, Often=1.0, Very often=1.0]
    - How often do you have difficulty getting things in order when you have to do a task that requires organization?  — [Never=0.0, Rarely=0.0, Sometimes=1.0, Often=1.0, Very often=1.0]
    - How often do you have problems remembering appointments or obligations?  — [Never=0.0, Rarely=0.0, Sometimes=1.0, Often=1.0, Very often=1.0]
    - When you have a task that requires a lot of thought, how often do you avoid or delay getting started?  — [Never=0.0, Rarely=0.0, Sometimes=0.0, Often=1.0, Very often=1.0]
    - How often do you fidget or squirm with your hands or feet when you have to sit down for a long time?  — [Never=0.0, Rarely=0.0, Sometimes=0.0, Often=1.0, Very often=1.0]
    - How often do you feel overly active and compelled to do things, like you were driven by a motor?  — [Never=0.0, Rarely=0.0, Sometimes=0.0, Often=1.0, Very often=1.0]

### UCLA / ULS-8 — Loneliness

- **Items:** 8
- **Per-item scoring:** 2 answer scales across items (shown per item below)
- **Total score:** prorated sum of 8 items; 2 reverse-keyed
- **Auto-built:** yes (comp_ucla_loneliness)
- **Questions:**
    - I lack companionship  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - There is no one I can turn to  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - I am an outgoing person *(reverse-keyed)*  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0]
    - I feel left out  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - I feel isolated from others  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - I can find companionship when I want it *(reverse-keyed)*  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0]
    - I am unhappy being so withdrawn  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]
    - People are around me but not with me  — [Often=3.0, Sometime=2.0, Rarely=1.0, Never=0.0, Sometimes=2.0]

### Everyday Discrimination Scale

- **Items:** 9
- **Per-item scoring:** Almost everyday = 6, At least once a week = 5, A few times a month = 4, A few times a year = 3, Less than once a year = 2, Never = 1
- **Total score:** prorated sum of 9 items; no reverse-keyed items
- **Auto-built:** yes (comp_everyday_discrimination)
- **Questions:**
    - You are treated with less courtesy than other people are.
    - You are treated with less respect than other people are.
    - You receive poorer service than other people at restaurants or stores.
    - People act as if they are afraid of you.
    - People act as if they're better than you are.
    - You are called names or insulted.
    - You are threatened or harassed.
    - People act as if they think you are not smart.
    - People act as if they think you are dishonest.

### MOS Social Support (RAND) + Tangible subscale

- **Items:** 9
- **Per-item scoring:** None of the time = 1, A little of the time = 2, Some of the time = 3, Most of the time = 4, All of the time = 5
- **Total score:** prorated sum of 9 items; no reverse-keyed items
- **Auto-built:** yes (comp_social_support)
- **Questions:**
    - Someone to help you if you were confined to bed
    - Someone to take you to the doctor if you needed it
    - Someone to prepare your meals if you were unable to do it yourself
    - Someone to help with daily chores if you were sick
    - Someone to take you to the doctor if you need it
    - Someone to have a good time with
    - Someone to turn to for suggestions about how to deal with a personal problem
    - Someone who understands your problems
    - Someone to love and make you feel wanted

### Neighborhood, walkability & food-insecurity composites

Built directly from the survey items (reusing their ordinal scores), because the scoring sheet groups these with mixed item valence. Opposite-valence items are reverse-keyed.

#### comp_social_cohesion

- Neighborhood social cohesion; higher = more cohesion.
- **Items:** 4; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - People around here are willing to help their neighbors.
    - People in my neighborhood generally get along with each other.
    - People in my neighborhood can be trusted.
    - People in my neighborhood share the same values.

#### comp_neighborhood_disorder

- Perceived neighborhood disorder (order items reversed); higher = more disorder.
- **Items:** 13; **reverse-keyed:** 4; prorated sum
- **Questions:**
    - There is a lot of graffiti in my neighborhood.
    - My neighborhood is noisy.
    - Vandalism is common in my neighborhood.
    - There are lot of abandoned buildings in my neighborhood.
    - There are too many people hanging around on the streets near my home.
    - There is a lot of crime in my neighborhood.
    - There is too much drug use in my neighborhood.
    - There is too much alcohol use in my neighborhood.
    - I'm always having trouble with my neighbors.
    - My neighborhood is clean. *(reverse-keyed)*
    - People in my neighborhood take good care of their houses and apartments. *(reverse-keyed)*
    - In my neighborhood, people watch out for each other. *(reverse-keyed)*
    - My neighborhood is safe. *(reverse-keyed)*

#### comp_neighborhood_physical_disorder

- Physical disorder subscale (order items reversed); higher = more disorder.
- **Items:** 6; **reverse-keyed:** 2; prorated sum
- **Questions:**
    - There is a lot of graffiti in my neighborhood.
    - My neighborhood is noisy.
    - Vandalism is common in my neighborhood.
    - There are lot of abandoned buildings in my neighborhood.
    - My neighborhood is clean. *(reverse-keyed)*
    - People in my neighborhood take good care of their houses and apartments. *(reverse-keyed)*

#### comp_neighborhood_social_disorder

- Social disorder subscale (order items reversed); higher = more disorder.
- **Items:** 7; **reverse-keyed:** 2; prorated sum
- **Questions:**
    - There are too many people hanging around on the streets near my home.
    - There is a lot of crime in my neighborhood.
    - There is too much drug use in my neighborhood.
    - There is too much alcohol use in my neighborhood.
    - I'm always having trouble with my neighbors.
    - In my neighborhood, people watch out for each other. *(reverse-keyed)*
    - My neighborhood is safe. *(reverse-keyed)*

#### comp_neighborhood_walkability

- PANES neighborhood walkability (crime-safety items reversed); higher = more walkable.
- **Items:** 7; **reverse-keyed:** 2; prorated sum
- **Questions:**
    - Many shops, stores, markets or other places to buy things I need are within easy walking distance of my home. Would you say that you...
    - It is within a 10-15 minute walk to a transit stop (such as bus, train, trolley, or tram) from my home. Would you say that you...
    - There are sidewalks on most of the streets in my neighborhood. Would you say that you...
    - There are facilities to bicycle in or near my neighborhood, such as special lanes, separate paths or trails, or shared use paths for cycles and pedestrians. Would you say that you...
    - My neighborhood has several free or low-cost recreation facilities, such as parks, walking trails, bike paths, recreation centers, playgrounds, public swimming pools, etc. Would you say that you...
    - The crime rate in my neighborhood makes it unsafe to go on walks at night. Would you say that you... *(reverse-keyed)*
    - The crime rate in my neighborhood makes it unsafe to go on walks during the day. Would you say that you... *(reverse-keyed)*

#### comp_hunger_vital_sign

- Hunger Vital Sign food-insecurity screener; higher = more food insecurity.
- **Items:** 2; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - Within the past 12 months, we worried whether our food would run out before we got money to buy more.
    - Within the past 12 months, the food we bought just didn't last and we didn't have money to get more.

#### comp_ptsd_pcl

- PTSD symptoms (abbreviated PCL, 5 items, 0-4 each); higher = more symptoms.
- **Items:** 5; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - In the past month, have you had repeated, disturbing memories, thoughts, or images of a stressful experience from the past?
    - In the past month, have you felt very upset when something reminded you of a stressful experience from the past?
    - In the past month, have you avoided activities or situations because they reminded you of a stressful experience from the past?
    - In the past month, have you felt distant or cut off from other people?
    - In the past month, have you felt irritable or had angry outbursts?

#### comp_subjective_wellbeing

- Subjective well-being (happiness + life meaning, UKB-style); higher = greater well-being.
- **Items:** 2; **reverse-keyed:** 0; prorated sum
- **Questions:**
    - In general, how happy are you?
    - To what extent do you feel your life to be meaningful?
