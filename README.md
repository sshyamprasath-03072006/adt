# 📚 THE ULTIMATE FOUNDATIONS OF DATA SCIENCE STUDY GUIDE
## 🎯 100/100 CENTUM GUARANTEE

---

## 📖 TABLE OF CONTENTS

**UNIT 1: INTRODUCTION AND DESCRIBING DATA** (Pages 1-50)
**UNIT 2: DESCRIBING RELATIONSHIPS AND PYTHON LIBRARIES** (Pages 51-120)
**UNIT 3: DATA WRANGLING AND DATA VISUALIZATION** (Pages 121-200)

---

# 🎓 UNIT 1: INTRODUCTION AND DESCRIBING DATA

## 1.1 DATA SCIENCE: FOUNDATIONS

### 1.1.1 Definition of Data Science

**Data Science** is the domain of study that deals with vast volumes of data using modern tools and techniques (like machine learning) to find unseen patterns, derive meaningful information, and make business decisions.

### 1.1.2 Benefits and Uses of Data Science

#### **Benefits:**

1. **Better Decision Making** — Data Science enables organizations to make data-driven decisions by analyzing patterns and trends hidden in large datasets, reducing guesswork.

2. **Identifying Target Audiences** — Businesses use data science to segment and identify their most valuable customer groups, enabling more effective marketing and resource allocation.

3. **Identifies Opportunities** — Helps discover new business opportunities through pattern recognition.

4. **Improves Performance** — Optimizes business processes and operational efficiency.

#### **Uses:**
- Predicts future outcomes
- Lowers risks (e.g., driverless cars)
- Customizes therapeutics (genomics)
- Fraud detection
- Recommendation systems

---

## 1.2 FACETS OF DATA

Data exists in multiple forms called **facets**. Understanding these is crucial for data science:

### 1. **Structured Data**
- Organized in rows and columns format
- Examples: Excel spreadsheets, SQL databases
- Easily searchable and understood by computers
- Makes data analysis straightforward

### 2. **Unstructured Data**
- Does not follow a specified format or structural rules
- Makes up **>80% of all data**
- Examples: Emails, PDFs, customer feedback, images, videos
- Requires specialized tools for analysis

### 3. **Semi-structured Data**
- Contains tags/markers to separate elements
- Examples: JSON, XML files
- Has some organizational properties but not rigid structure

### 4. **Natural Language**
- Human speech and text
- Requires **Natural Language Processing (NLP)** to understand
- Examples: Social media posts, emails, documents

### 5. **Machine-Generated Data**
- Created without human interaction by computer processes
- Examples: Web server logs, IoT sensor data, system logs
- Typically high-volume and time-stamped

### 6. **Graph-Based (Network) Data**
- Describes relationships (edges) between entities (nodes)
- Examples: Social network followers, organizational hierarchies
- Requires graph databases and algorithms

### 7. **Audio, Image, and Video**
- Multimedia data
- Trivial for humans to interpret
- Poses significant challenges for computer analysis
- Requires deep learning and computer vision techniques

---

## 1.3 THE DATA SCIENCE PROCESS

The Data Science Process is a structured lifecycle followed to convert raw data into meaningful insights and actionable decisions. It consists of **6 key steps**:

### **Step 1: Setting the Research Goal**

This is the **foundation** of any data science project.

**Key Activities:**
- Define the **what, how, and why** of the project
- Ensure all stakeholders are aligned
- Create a **Project Charter**

**Project Charter Contains:**
- Clear research goal and objectives
- Project mission and context
- Analysis plan and required resources
- Proof of concept and deliverables
- Realistic timeline

**⚠️ Warning:** Without a well-defined research goal, the project can drift in the wrong direction, wasting time and resources.

---

### **Step 2: Retrieving Data**

Once the goal is defined, gather the necessary data.

**Data Sources:**

**Internal Sources:**
- Company databases
- Data marts
- Data warehouses
- Data lakes (historical + current data stored internally)

**External Sources:**
- Purchased or licensed data from third-party providers
- Examples: Nielsen, Twitter APIs, government datasets

**Quality Check:** Data quality checks are performed at this stage to prevent easily avoidable errors from corrupting the rest of the project.

---

### **Step 3: Data Preparation (CRUCIAL - 80% of Time!)**

This is the **most time-consuming step** — data scientists typically spend **up to 80%** of the project's time here.

#### **A. Data Cleansing:**

1. **Fixing errors and typos**
   - Correct spelling mistakes
   - Fix formatting inconsistencies

2. **Handling Missing Values:**
   - **Imputation**: Fill in missing values
     - Fill with mean/median/mode
     - Forward fill (ffill): Use previous value
     - Backward fill (bfill): Use next value
   - **Deletion**: Remove rows/columns with missing data

3. **Handling Outliers:**
   - An extreme data point that differs significantly from other observations
   - **Options:**
     - Cap them (set maximum threshold)
     - Remove them entirely
     - Transform the data

4. **Removing Duplicates:**
   - Identify and delete duplicate records
   - Prevents bias in analysis

#### **B. Data Transformation:**

1. **Aggregating data** to a suitable level
2. **Deriving new computed measures** (e.g., profit margin from sales and cost)
3. **Creating dummy variables** for categorical data
4. **Joining multiple datasets** into one unified table
5. **Converting formats** (e.g., text to numbers)

---

### **Step 4: Data Exploration (EDA - Exploratory Data Analysis)**

**EDA** means diving deeper into the cleaned data using descriptive statistics and visualization techniques to discover hidden patterns, relationships, and anomalies before building models.

**Techniques Used:**

1. **Simple Graphs:**
   - Histograms (distribution of single variable)
   - Density plots
   - Box plots

2. **Combined Graphs:**
   - Scatter plots (relationships between two variables)

3. **Brushing and Linking:**
   - Selecting observations in one plot automatically highlights those same observations in all other plots simultaneously
   - Very useful for spotting multivariate patterns

**Purpose:** EDA helps the analyst understand the data's nature before committing to a specific model.

---

### **Step 5: Data Modeling**

This is the **core analytical step** where machine learning and statistical techniques are applied to achieve the research goal.

**Key Sub-tasks:**

1. **Model Selection:**
   - Choose the algorithm based on:
     - Performance requirements
     - Ease of implementation
     - Explainability
     - Maintenance cost
   - Examples: Linear regression, decision trees, clustering

2. **Execution:**
   - Writing code using Python libraries:
     - Scikit-learn
     - StatsModels
     - NumPy, Pandas

3. **Model Evaluation:**
   - Use a **Holdout Sample strategy**
   - Typical split: **80% training / 20% testing**
   - Train model on training data
   - Evaluate on unseen test data
   - Error measures:
     - Mean Square Error (MSE)
     - R² (coefficient of determination)
   - Ensures the model generalizes well

---

### **Step 6: Presentation and Automation**

The final step involves communicating findings to non-technical stakeholders and making the model operational.

**A. Presentation:**
- Create clear visual reports
- Build dashboards
- Prepare presentations
- Translate complex findings into actionable business recommendations

**B. Automation:**
- Set up automated pipelines
- Model can score new data automatically
- Update reports without manual intervention
- Refresh Excel/PowerPoint outputs automatically

---

### **Summary Table: Data Science Process**

| Step | Key Activity | Output |
|------|-------------|--------|
| 1. Research Goal | Define the problem | Project Charter |
| 2. Data Retrieval | Collect internal/external data | Raw dataset |
| 3. Data Preparation | Cleanse and transform | Clean dataset |
| 4. EDA | Visual & statistical exploration | Insights, patterns |
| 5. Data Modeling | Build, train, and test model | Predictive model |
| 6. Presentation | Communicate & automate | Reports, dashboards |

---

## 1.4 DATA MINING vs DATA WAREHOUSING

### **DATA WAREHOUSING**

**Definition:**
A data warehouse is a **subject-oriented, integrated, time-varying, and non-volatile** database system designed specifically for analytical analysis rather than transactional (day-to-day) operations.

**Key Characteristics:**

1. **Subject-Oriented:**
   - Organized around key business subjects like sales, customers, or products
   - NOT organized around operational processes

2. **Integrated:**
   - Combines data from multiple heterogeneous sources
   - Converts into a consistent unified format

3. **Time-Varying:**
   - Stores historical data over time
   - Enables trend analysis

4. **Non-Volatile:**
   - Data is loaded and read
   - NOT frequently updated or deleted

**Process - ETL (Extract, Transform, Load):**

```
┌─────────────┐    ┌──────────────┐    ┌─────────────┐
│   EXTRACT   │ -> │  TRANSFORM   │ -> │    LOAD     │
└─────────────┘    └──────────────┘    └─────────────┘
      │                    │                    │
      v                    v                    v
Multiple sources    Clean, standardize,   Warehouse for
(Sales, CRM, etc.)  restructure data      querying
```

1. **Extract:** Pull data from multiple source systems
2. **Transform:** Clean, standardize, and restructure the data
3. **Load:** Store the processed data into the warehouse for querying

**Advantages:**
- Easy to understand and query
- Provides continuous updates
- Stores historical data for trend evaluation
- Enables faster, organization-wide reporting

**Disadvantages:**
- Risk of accumulating irrelevant or outdated data
- Cleansing data arriving from multiple heterogeneous sources is complex

**Managed by:** Data Engineers / IT teams

**Example:** Consolidating sales from 50 stores into one server

---

### **DATA MINING**

**Definition:**
Data Mining is the process of analyzing large datasets using **AI, machine learning, and statistical techniques** to discover meaningful patterns, relationships, and trends, and to predict future outcomes.

**Key Techniques:**

1. **Classification:**
   - Assigning items into predefined categories
   - Example: Spam vs. not spam email detection

2. **Clustering:**
   - Grouping similar data points without predefined labels
   - Example: Customer segmentation

3. **Association Rule Mining:**
   - Finding "if-then" patterns
   - Example: Market basket analysis
   - "Customers who buy bread also buy butter"

4. **Prediction/Regression:**
   - Forecasting continuous numerical outcomes
   - Example: Predicting house prices

**Advantages:**
- Enables fault detection and fraud prevention
- Cost-effective once implemented
- Produces easily accessible, actionable knowledge

**Disadvantages:**
- Not 100% accurate; incorrect patterns can cause breaches
- Resource-heavy: training and implementing models requires significant computation

**Managed by:** Business users and Data Scientists

**Example:** Finding "customers who buy X also buy Y"

---

### **COMPARISON TABLE: Data Warehousing vs Data Mining**

| Feature | Data Warehousing | Data Mining |
|---------|-----------------|-------------|
| **Definition** | Database system for analytical analysis | Process of analyzing data to find patterns |
| **Process** | ETL (Extract, Transform, Load) | AI, ML, Statistics |
| **Purpose** | Store & report on historical data | Discover hidden patterns & predict outcomes |
| **Data Handling** | Pools all relevant data together | Extracts knowledge from large datasets |
| **Managed By** | Data Engineers | Data Scientists / Business Users |
| **Functionality** | Subject-oriented, integrated, time-varying, non-volatile | AI, statistics, databases, machine learning |
| **Task** | Extracting and storing data for efficient reporting | Pattern recognition logic to find patterns |
| **Uses** | Makes reporting easier and faster | Aids in identification of access patterns |
| **Example** | CRM system integration | Customer purchasing behavior analysis |

---

## 1.5 TYPES OF DATA

### **1. Qualitative vs Quantitative Data**

**Qualitative (Categorical) Data:**
- Describes qualities or characteristics
- Non-numeric
- Examples: Gender, color, brand name, yes/no responses

**Quantitative (Numerical) Data:**
- Represents measurable quantities
- Can be counted or measured
- Two types:

  **a) Discrete Data:**
  - Countable, whole numbers only
  - Cannot have decimals
  - Examples: Number of students, number of cars, number of rooms
  
  **b) Continuous Data:**
  - Measurable, infinite decimals allowed
  - Can take any value within a range
  - Examples: Height, temperature, weight, time

---

### **2. Four Scales of Measurement (NOIR)**

**N - Nominal Scale:**
- Categories with **NO order**
- Just labels or names
- Examples: Gender (Male/Female), Colors (Red/Blue/Green), Country names
- Cannot perform mathematical operations

**O - Ordinal Scale:**
- Categories **WITH order**
- Exact gaps between them are **unknown**
- Examples: 
  - 1-5 star ratings (⭐⭐⭐⭐⭐)
  - Education level (High School < Bachelor's < Master's < PhD)
  - Satisfaction (Very Unsatisfied < Unsatisfied < Neutral < Satisfied < Very Satisfied)

**I - Interval Scale:**
- Ordered with **equal intervals**
- **NO true zero** point
- Zero does NOT mean "absence of value"
- Examples: 
  - Temperature in °C (0°C does not mean "no temperature")
  - Calendar years
- Can add/subtract but not multiply/divide meaningfully

**R - Ratio Scale:**
- Ordered with equal intervals
- **WITH a true zero point**
- Zero means complete absence of the quantity
- Examples: 
  - Weight (0 kg means no weight)
  - Income (₹0 means no income)
  - Age, Height, Distance
- All mathematical operations valid

---

## 1.6 TYPES OF VARIABLES

### **1. Independent Variable (IV)**
- The hypothesized **cause**
- The variable you **manipulate** or use to **predict**
- Also called: Predictor variable, Explanatory variable
- Plotted on **X-axis** (horizontal)
- Example: Study hours (used to predict exam scores)

### **2. Dependent Variable (DV)**
- The **effect** or **outcome**
- The variable being **measured**
- Also called: Response variable, Outcome variable
- Plotted on **Y-axis** (vertical)
- Example: Exam scores (depends on study hours)

### **3. Control Variable**
- Held **constant** to prevent bias
- Ensures fair comparison
- Example: Keeping room temperature the same during an experiment

### **4. Confounding Variable**
- A **hidden variable** that affects both the IV and DV
- Causes a **false relationship**
- Must be identified and controlled
- Example: In studying relationship between ice cream sales and drowning deaths, **temperature** is a confounding variable (hot weather increases both)

---

## 1.7 DESCRIPTIVE STATISTICS

### **1.7.1 Measures of Central Tendency**

These describe where the "center" of the data lies.

#### **A. Mean (Arithmetic Average)**

**Formula:**
```
        Σx
μ = ───────
        n
```

Where:
- μ = Mean
- Σx = Sum of all values
- n = Number of values

**Characteristics:**
- Best for **symmetric data**
- **Sensitive to outliers** (extreme values can skew it)
- Uses all data points

**Example 1:** Calculate mean for data: {40, 45, 50, 55, 60}

**Step-by-Step Solution:**
```
Step 1: Add all values
Sum = 40 + 45 + 50 + 55 + 60 = 250

Step 2: Count number of values
n = 5

Step 3: Apply formula
      250
μ = ───── = 50
       5

Answer: Mean = 50
```

**Example 2:** Data Science Application
```
Monthly sales (₹ lakhs): [120, 75, 180, 145, 60]

Sum = 120 + 75 + 180 + 145 + 60 = 580
n = 5

Mean sales = 580/5 = ₹116 lakhs
```

---

#### **B. Median**

**Definition:** The exact **middle value** when data is ordered from smallest to largest.

**How to Find:**
1. Arrange data in **ascending order**
2. If **n is odd**: Median is the middle value
3. If **n is even**: Median is the average of two middle values

**Characteristics:**
- **NOT affected by outliers**
- Preferred for **skewed data**
- More robust than mean

**Example 1:** Odd number of values
```
Data: {40, 45, 50, 55, 60}

Step 1: Already ordered
Step 2: n = 5 (odd)
Step 3: Middle position = (5+1)/2 = 3rd value
Step 4: Median = 50

Answer: Median = 50
```

**Example 2:** Even number of values
```
Data: {10, 20, 30, 40}

Step 1: Already ordered
Step 2: n = 4 (even)
Step 3: Two middle positions = 2nd and 3rd values
Step 4: Median = (20 + 30)/2 = 25

Answer: Median = 25
```

**Example 3:** With Outliers
```
Salaries (₹): [30000, 32000, 35000, 38000, 200000]

Mean = (30000+32000+35000+38000+200000)/5 = ₹67,000
Median = 35000 (3rd value)

The outlier (₹200,000) heavily skews the mean,
but median remains representative!
```

---

#### **C. Mode**

**Definition:** The **most frequently occurring value** in a dataset.

**Characteristics:**
- Best for **nominal/categorical data**
- Dataset can be:
  - **Unimodal**: One mode
  - **Bimodal**: Two modes
  - **Multimodal**: More than two modes
  - **No mode**: All values occur equally

**Example 1:** Single mode
```
Data: {10, 20, 20, 30, 40}

Frequency:
10 → 1 time
20 → 2 times ← Most frequent
30 → 1 time
40 → 1 time

Mode = 20
```

**Example 2:** Bimodal
```
Data: {5, 5, 10, 10, 15, 20}

Frequency:
5 → 2 times
10 → 2 times
15 → 1 time
20 → 1 time

Modes = 5 and 10 (Bimodal)
```

**Example 3:** No mode
```
Data: {1, 2, 3, 4, 5}

All values occur once → No mode
```

---

### **1.7.2 Measures of Variability (Dispersion)**

These describe how **spread out** the data is.

#### **A. Range**

**Formula:**
```
Range = Maximum value - Minimum value
```

**Characteristics:**
- Simplest measure
- **Heavily affected by outliers**
- Only uses two values

**Example:**
```
Data: {40, 45, 50, 55, 60}

Maximum = 60
Minimum = 40

Range = 60 - 40 = 20
```

---

#### **B. Variance (σ²)**

**Definition:** Average of **squared deviations** from the mean.

**Formula (Population):**
```
           Σ(x - μ)²
σ² = ───────────────
              n
```

Where:
- σ² = Variance
- x = Each value
- μ = Mean
- n = Number of values

**Characteristics:**
- Uses all data points
- Squared units (hard to interpret directly)
- Always positive

**Example 1:** Calculate variance for data: {40, 45, 50, 55, 60}

**Step-by-Step Solution:**
```
Step 1: Calculate mean
μ = (40+45+50+55+60)/5 = 250/5 = 50

Step 2: Find deviation of each value from mean
(40 - 50) = -10
(45 - 50) = -5
(50 - 50) = 0
(55 - 50) = 5
(60 - 50) = 10

Step 3: Square each deviation
(-10)² = 100
(-5)²  = 25
(0)²   = 0
(5)²   = 25
(10)²  = 100

Step 4: Sum of squared deviations
Σ(x - μ)² = 100 + 25 + 0 + 25 + 100 = 250

Step 5: Divide by n
       250
σ² = ───── = 50
        5

Answer: Variance = 50
```

---

#### **C. Standard Deviation (σ)**

**Definition:** Square root of variance. Shows the average spread of data around the mean in **original units**.

**Formula:**
```
σ = √σ² = √[Σ(x - μ)² / n]
```

**Characteristics:**
- Expressed in same units as data
- Most widely used measure of spread
- Smaller σ = data clustered around mean
- Larger σ = data more spread out

**Example:** Using previous variance calculation
```
Data: {40, 45, 50, 55, 60}
Variance (σ²) = 50

Standard Deviation:
σ = √50 = 7.07

Interpretation: On average, values deviate 
from the mean by about 7.07 units.
```

---

#### **D. Coefficient of Variation (CV)**

**Formula:**
```
        σ
CV = ───── × 100%
        μ
```

**Characteristics:**
- **Relative measure** of variability
- Allows comparison across datasets with **different units**
- Unitless (percentage)

**Interpretation:**
- **CV < 10%** → Low variability (consistent data)
- **CV > 30%** → High variability (highly variable data)

**Example 1:** Comparing two datasets
```
Dataset A: Mean = 100, SD = 10
Dataset B: Mean = 50, SD = 10

CV_A = (10/100) × 100% = 10%
CV_B = (10/50) × 100% = 20%

Even though both have same SD, 
Dataset B is MORE variable relative to its mean!
```

**Example 2:** Real-world application
```
Stock A: Average return = 12%, SD = 8%
Stock B: Average return = 8%, SD = 5%

CV_A = (8/12) × 100% = 66.7%
CV_B = (5/8) × 100% = 62.5%

Stock A has slightly higher relative risk.
```

---

### **1.7.3 Important Theory: Linear Transformation**

**Concept:** If you add or multiply a constant to every value in a dataset, how do statistics change?

**Case 1: Adding a Constant**
```
If you add a constant (e.g., +5) to every value:

✓ Mean increases by that constant
✓ Median increases by that constant
✗ Variance remains SAME
✗ Standard Deviation remains SAME

Why? The spread/distance between numbers hasn't changed,
only the location shifted.
```

**Example:**
```
Original data: {10, 20, 30}
Mean = 20, Variance = 66.67, SD = 8.16

Add 5 to each: {15, 25, 35}
New Mean = 25 (increased by 5)
New Variance = 66.67 (SAME)
New SD = 8.16 (SAME)
```

**Case 2: Multiplying by a Constant**
```
If you multiply every value by a constant (e.g., ×2):

✓ Mean is multiplied by that constant
✓ Median is multiplied by that constant
✓ Variance is multiplied by constant²
✓ Standard Deviation is multiplied by constant
```

---

### **1.7.4 Distribution Shapes**

#### **A. Skewness (Symmetry)**

Measures the **lack of symmetry** in a distribution.

**1. Symmetric Distribution:**
```
           ___
         _/   \_
        /       \
       /         \
      /           \
     /             \
    /_______________\
    
Mean = Median = Mode
```

**Characteristics:**
- Perfectly balanced
- Bell-shaped (Normal distribution)
- Example: Heights of adults

**2. Right Skewed (Positive Skew):**
```
           _
         _/ \_
        /     \__
       /         \__
      /             \___
     /                  \__
    /______________________\___
    
Mean > Median > Mode
(Long tail on right)
```

**Characteristics:**
- Long tail on the **right** side
- Mean pulled toward tail
- Example: Income distribution (few very high earners)

**3. Left Skewed (Negative Skew):**
```
   _
 _/ \_
/     \__
         \__
            \___
                \__
___________________\___

Mean < Median < Mode
(Long tail on left)
```

**Characteristics:**
- Long tail on the **left** side
- Mean pulled toward tail
- Example: Age at retirement (most retire late, few retire very early)

---

#### **B. Kurtosis (Peakedness)**

Measures the **"peakedness"** or "tailedness" of a distribution.

**1. Mesokurtic (Normal):**
```
           ___
         _/   \_
        /       \
       /         \
      /           \
     /             \
    /_______________\
    
Normal bell shape
Kurtosis = 3
```

**2. Leptokurtic (High Peak):**
```
            _
           / \
          /   \
         /     \
        /       \
       /         \
      /           \
     /             \
    /_______________\
    
More peaked, heavier tails
Kurtosis > 3
More outliers than normal
```

**3. Platykurtic (Flat):**
```
         _______
        /       \
       /         \
      /           \
     /             \
    /_______________\
    
Flatter, lighter tails
Kurtosis < 3
Fewer outliers than normal
```

---

## 1.8 VISUALIZING DATA

### **1.8.1 Frequency Distribution Table**

A systematic way to show the **frequency (count)** of values.

**Types:**

**A. Ungrouped Frequency Distribution:**
Shows individual numbers.

**Example:**
```
Marks    Frequency
  10         2
  15         5
  20         8
  25         3
  30         2
```

**B. Grouped Frequency Distribution:**
Uses class intervals.

**Example:**
```
Class Interval    Frequency
    10-14             3
    15-19             7
    20-24            12
    25-29             5
    30-34             3
```

---

### **Example: Constructing a Grouped Frequency Distribution**

**Problem:** Given 20 student marks (range 10 to 25). Create class intervals, frequency table, and find mean.

**Data:** 12, 15, 18, 14, 20, 22, 11, 16, 19, 23, 13, 17, 21, 10, 24, 15, 18, 20, 16, 19

**Step-by-Step Solution:**

```
Step 1: Determine Class Intervals
Range = Max - Min = 24 - 10 = 14
Let's use 4 class intervals of width 4

Class Intervals: 10-13, 14-17, 18-21, 22-25

Step 2: Create Tally and Count Frequency

Class      Tally          Frequency (f)
10-13     |||| |              5
14-17     |||| ||             7
18-21     |||| ||             6
22-25     ||                  2
                     ─────────────
Total                       20

Step 3: Find Midpoint (x) of Each Interval

Class      f      x      f × x
10-13      5     11.5     57.5
14-17      7     15.5    108.5
18-21      6     19.5    117.0
22-25      2     23.5     47.0
                 ─────────────
                Σf=20   Σ(f×x)=330

Step 4: Calculate Grouped Mean

         Σ(f × x)    330
x̄ = ──────────── = ──── = 16.5
           Σf        20

Answer: Mean = 16.5
```

---

### **1.8.2 Histogram**

**Definition:** Contiguous (touching) bars representing **continuous data** frequencies.

**Characteristics:**
- Bars **touch** each other (no gaps)
- X-axis: Class intervals
- Y-axis: Frequency
- Area of bar represents frequency

**ASCII Diagram:**
```
Frequency
   12 |           ████
   10 |           ████
    8 |     ████  ████
    6 |     ████  ████  ████
    4 |     ████  ████  ████
    2 | ████ ████  ████  ████
    0 |███████████████████████
      10   15   20   25   30  Marks
```

**When to Use:**
- Checking normality
- Spotting skewness
- Identifying outliers
- Understanding distribution shape

---

### **1.8.3 Frequency Polygon**

**Definition:** Line graph connecting the **midpoints** of histogram bars.

**Characteristics:**
- Good for **comparing multiple distributions**
- Smoother appearance than histogram
- Can overlay multiple polygons

**ASCII Diagram:**
```
Frequency
   12 |           ●
   10 |           ●
    8 |     ●     ●
    6 |     ●     ●     ●
    4 |     ●     ●     ●
    2 | ●   ●     ●     ●
    0 |●____●_____●_____●____●
      10   15   20   25   30
```

---

### **1.8.4 Bar Chart**

**Definition:** Bars that **do NOT touch**. Used for **categorical/nominal** data.

**Characteristics:**
- Gaps between bars
- Each bar represents a category
- Height/length represents value

**ASCII Diagram:**
```
Sales (₹K)
  100 |  ████
   80 |  ████        ████
   60 |  ████  ████  ████
   40 |  ████  ████  ████  ████
   20 |  ████  ████  ████  ████
    0 |__████__████__████__████__
        North South East  West
         Region
```

**When to Use:**
- Comparing sales across regions
- Comparing scores across subjects
- Comparing categories

---

### **1.8.5 Box Plot (Box-and-Whisker Plot)**

**Definition:** Visualizes **quartiles (Q1, Median, Q3)** and outliers.

**Components:**
- **Q1 (First Quartile):** 25th percentile
- **Median (Q2):** 50th percentile
- **Q3 (Third Quartile):** 75th percentile
- **IQR (Interquartile Range):** Q3 - Q1
- **Whiskers:** Extend to minimum and maximum (excluding outliers)
- **Outliers:** Points beyond whiskers

**ASCII Diagram:**
```
    Outlier
       ○
       │
       │
  Min ─┼───
       │   │
  Q1 ──┼───┤
       │   │
Median─┼───┤
       │   │
  Q3 ──┼───┤
       │   │
  Max ─┼───
       │
       ○
    Outlier
```

**When to Use:**
- Comparing distributions across groups
- Identifying outliers
- Showing spread and central tendency

---

### **1.8.6 Pie Chart**

**Definition:** Shows **proportions** as slices of a circle.

**Characteristics:**
- Each slice represents a category
- Size proportional to percentage
- Total = 100%

**ASCII Diagram:**
```
         Sales Distribution
         
           _______
        .-'       '-.
      .'   North    '.
     /     40%        \
    |                  |
    |   West    East   |
    |   15%     25%    |
     \                 /
      '.   South     .'
        '-. 20%   .-'
           '-----'
```

**When to Use:**
- Showing market share
- Displaying budget allocation
- Percentage breakdown

---

### **1.8.7 Line Graph**

**Definition:** Shows **trends over time**.

**Characteristics:**
- X-axis: Time periods
- Y-axis: Measured values
- Points connected by lines

**ASCII Diagram:**
```
Sales (₹L)
  300 |                    ●
      |                 ╱
  250 |              ●
      |           ╱
  200 |        ●
      |     ╱
  150 |  ●
      |╱
  100 |
      └─────────────────────────
      2019 2020 2021 2022 2023
           Year
```

**When to Use:**
- Revenue growth over years
- Stock price trends
- Temperature changes over months

---

## 1.9 NORMAL DISTRIBUTION & Z-SCORES

### **1.9.1 Normal Distribution**

**Definition:** A **symmetric, bell-shaped** distribution where:
- Mean = Median = Mode
- Tails never touch the X-axis (theoretically infinite)
- Completely defined by mean (μ) and standard deviation (σ)

**ASCII Diagram:**
```
            ___
          _/   \_
         /       \
        /         \
       /           \
      /             \
     /               \
    /                 \
   /___________________\
   
 μ-3σ  μ-2σ  μ-σ   μ   μ+σ  μ+2σ  μ+3σ
```

**Characteristics:**
- Symmetric around mean
- 68-95-99.7 rule applies
- Most common distribution in nature
- Foundation for many statistical tests

---

### **1.9.2 Empirical Rule (68-95-99.7 Rule)**

For a **normal distribution**:

```
            68%
          ←─────→
            ___
          _/   \_
         /       \
        /         \
       /           \
      /             \
     /               \
    /                 \
   /___________________\
   
 μ-3σ  μ-2σ  μ-σ   μ   μ+σ  μ+2σ  μ+3σ
 ←───────── 95% ─────────→
 ←─────────────── 99.7% ───────────────→
```

**The Rule States:**

1. **68%** of data falls within **±1 Standard Deviation** of the mean
   - Range: (μ - 1σ) to (μ + 1σ)

2. **95%** of data falls within **±2 Standard Deviations** of the mean
   - Range: (μ - 2σ) to (μ + 2σ)

3. **99.7%** of data falls within **±3 Standard Deviations** of the mean
   - Range: (μ - 3σ) to (μ + 3σ)

---

### **Example Problems: Empirical Rule**

**Example 1:** Basic Application

**Problem:** Dataset has Mean = 50, SD = 8.
a) What % of data lies between 42 and 58?
b) What % lies beyond 66?

**Step-by-Step Solution:**

```
Given:
μ = 50
σ = 8

Part (a): % between 42 and 58

Step 1: Convert to standard deviations
42 = 50 - 8 = μ - 1σ
58 = 50 + 8 = μ + 1σ

Step 2: Apply empirical rule
Range is μ ± 1σ
Therefore, 68% of data lies here

Answer (a): 68%

Part (b): % beyond 66

Step 1: Convert to standard deviations
66 = 50 + 16 = 50 + 2(8) = μ + 2σ

Step 2: Apply empirical rule
95% lies between μ - 2σ and μ + 2σ
That is, between 34 and 66

Step 3: Calculate remaining percentage
100% - 95% = 5% lies outside this range

Step 4: Divide by 2 (symmetry)
Since distribution is symmetric:
5%/2 = 2.5% lies below 34
2.5% lies above 66

Answer (b): 2.5% lies beyond 66
```

---

**Example 2:** Advanced Application

**Problem:** IQ scores are normally distributed with μ = 100, σ = 15.
What percentage of people have IQ:
a) Between 85 and 115?
b) Above 130?
c) Between 70 and 130?

**Solution:**

```
Given: μ = 100, σ = 15

Part (a): Between 85 and 115
85 = 100 - 15 = μ - 1σ
115 = 100 + 15 = μ + 1σ
Range = μ ± 1σ
Answer: 68%

Part (b): Above 130
130 = 100 + 30 = 100 + 2(15) = μ + 2σ

95% lies between μ - 2σ and μ + 2σ (70 to 130)
Remaining 5% is split equally in both tails
5%/2 = 2.5%

Answer: 2.5% have IQ above 130

Part (c): Between 70 and 130
70 = 100 - 30 = μ - 2σ
130 = 100 + 30 = μ + 2σ
Range = μ ± 2σ
Answer: 95%
```

---

### **1.9.3 Z-Score (Standard Score)**

**Definition:** A standardized score showing exactly **how many standard deviations** a value is from the mean.

**Formula:**
```
        X - μ
Z = ───────────
          σ
```

Where:
- Z = Z-score
- X = Raw score (observed value)
- μ = Mean
- σ = Standard deviation

**Interpretation:**
- **Z = 0**: Value equals the mean
- **Z > 0**: Value is above the mean
- **Z < 0**: Value is below the mean
- **Z = +2**: Value is 2 SDs above mean
- **Z = -1.5**: Value is 1.5 SDs below mean

**Advantages:**
- Allows comparison across different scales
- Identifies outliers (|Z| > 3)
- Used in hypothesis testing
- Converts any normal distribution to standard normal (μ=0, σ=1)

---

### **Example Problems: Z-Score**

**Example 1:** Basic Z-Score Calculation

**Problem:** Class mean = 60, SD = 5. Student scored 70. Find Z-score and interpret.

**Step-by-Step Solution:**

```
Given:
μ = 60
σ = 5
X = 70

Step 1: Apply formula
        X - μ      70 - 60      10
Z = ─────────── = ───────── = ──── = 2
          σ           5         5

Step 2: Interpretation
Z = +2.0 means:
- Student scored exactly 2 standard deviations 
  ABOVE the class mean
- This indicates excellent performance
- Approximately at 97.7th percentile
  (using standard normal table)

Answer: Z = 2.0 (Excellent performance)
```

---

**Example 2:** Comparing Performance Across Tests

**Problem:** 
- Math test: Score = 85, Class mean = 75, SD = 10
- English test: Score = 90, Class mean = 85, SD = 5

In which subject did the student perform better relative to classmates?

**Solution:**

```
Math:
        85 - 75      10
Z_M = ───────── = ──── = 1.0
          10        10

English:
        90 - 85       5
Z_E = ───────── = ──── = 1.0
           5         5

Interpretation:
Both Z-scores = 1.0
Student performed equally well relative to 
classmates in BOTH subjects!

Even though raw score was higher in English (90 vs 85),
the relative performance is the same.
```

---

**Example 3:** Identifying Outliers

**Problem:** Heights of students: μ = 165 cm, σ = 10 cm.
Is a height of 195 cm unusual?

**Solution:**

```
        195 - 165      30
Z = ─────────── = ──── = 3.0
          10          10

Interpretation:
Z = 3.0 means this height is 3 SDs above mean.
According to empirical rule, only 0.3% of data
lies beyond ±3σ.

Conclusion: YES, this is an outlier/unusual value.
Only about 0.15% of students would be this tall.
```

---

**Example 4:** Finding Raw Score from Z-Score

**Problem:** Test has μ = 70, σ = 8. What score corresponds to Z = 1.5?

**Solution:**

```
Formula rearrangement:
X = μ + (Z × σ)

X = 70 + (1.5 × 8)
X = 70 + 12
X = 82

Answer: A Z-score of 1.5 corresponds to a raw score of 82.
```

---

## 1.10 DESCRIBING DATA WITH TABLES AND GRAPHS - COMPREHENSIVE EXAMPLES

### **Example: Complete Data Analysis**

**Scenario:** A teacher has exam scores for 25 students:

**Raw Data:**
45, 52, 58, 61, 63, 65, 67, 68, 70, 72, 73, 75, 76, 78, 79, 80, 82, 84, 85, 87, 88, 90, 92, 95, 98

**Complete Analysis:**

```
STEP 1: Frequency Distribution

Class Interval    Frequency    Midpoint
40-49                1           44.5
50-59                2           54.5
60-69                5           64.5
70-79                7           74.5
80-89                6           84.5
90-99                4           94.5
                  ─────────
Total               25

STEP 2: Calculate Statistics

Mean:
Sum = 45+52+58+...+95+98 = 1893
n = 25
Mean = 1893/25 = 75.72

Median:
n = 25 (odd)
Median position = (25+1)/2 = 13th value
Median = 76

Mode:
All values appear once → No mode

Range:
Range = 98 - 45 = 53

Variance and Standard Deviation:
Σ(x - μ)² = 5847.84
Variance = 5847.84/25 = 233.91
SD = √233.91 = 15.29

Coefficient of Variation:
CV = (15.29/75.72) × 100% = 20.2%
Interpretation: Moderate variability

STEP 3: Quartiles

Q1 (25th percentile):
Position = 0.25 × 25 = 6.25 → 7th value
Q1 = 67

Q3 (75th percentile):
Position = 0.75 × 25 = 18.75 → 19th value
Q3 = 85

IQR = Q3 - Q1 = 85 - 67 = 18

STEP 4: Check for Outliers

Lower fence = Q1 - 1.5×IQR = 67 - 27 = 40
Upper fence = Q3 + 1.5×IQR = 85 + 27 = 112

All values between 40 and 112 → No outliers

STEP 5: Distribution Shape

Mean (75.72) ≈ Median (76)
→ Approximately symmetric distribution

STEP 6: Histogram Visualization

Frequency
  8 |           ████
  7 |           ████
  6 |                 ████
  5 |     ████
  4 |                       ████
  3 |
  2 |     ████
  1 | ████
  0 |███████████████████████████
    40  50  60  70  80  90  100
           Score Range

STEP 7: Box Plot

    ┌─────────┐
45──┤         ├──98
    └─────────┘
    67   76   85
    Q1  Med  Q3

Symmetric box, no outliers

CONCLUSION:
- Class performed moderately well (mean = 75.72)
- Scores are fairly consistent (CV = 20.2%)
- Distribution is approximately normal
- No extreme outliers
- Middle 50% of students scored between 67 and 85
```

---

## 📝 UNIT 1 SUMMARY CHECKLIST

✅ **Data Science Fundamentals**
- Definition and benefits understood
- 7 facets of data memorized
- 6-step data science process mastered

✅ **Data Mining vs Warehousing**
- Clear distinction made
- ETL process understood
- Use cases identified

✅ **Types of Data & Variables**
- Qualitative vs Quantitative
- NOIR scales mastered
- IV, DV, Control, Confounding variables

✅ **Descriptive Statistics**
- Mean, Median, Mode calculated
- Variance, SD, CV computed
- Linear transformation effects known

✅ **Distribution Shapes**
- Skewness identified
- Kurtosis types recognized
- Symmetry assessed

✅ **Data Visualization**
- Histogram, Bar chart, Pie chart
- Box plot, Line graph, Frequency polygon
- Appropriate use cases

✅ **Normal Distribution**
- Empirical rule (68-95-99.7)
- Z-score calculation
- Percentile interpretation

---

**END OF UNIT 1**

---

# 🎓 UNIT 2: DESCRIBING RELATIONSHIPS AND PYTHON LIBRARIES

## 2.1 CORRELATION

### **2.1.1 Definition**

**Correlation** is a statistical measure that describes the **strength and direction** of the **linear relationship** between two quantitative variables.

**Key Points:**
- Always **symmetric**: Correlation of X with Y = Correlation of Y with X
- **Does NOT imply causation**
- Range: **-1.0 to +1.0**
- Only measures **linear** relationships

---

### **2.1.2 Types of Correlation**

#### **1. Positive Correlation (r > 0)**

Both variables **increase together**.

**Strong Positive (r = 0.8 to 1.0):**
```
Y
│         ●
│       ●
│     ●
│   ●
│ ●
└───────────── X

Points tightly clustered in upward direction
```

**Weak Positive (r = 0.3 to 0.7):**
```
Y
│    ●   ●
│  ●  ●
│●    ●  ●
│  ●
└───────────── X

Points loosely scattered in upward trend
```

**Example:** Study hours and exam scores

---

#### **2. Negative Correlation (r < 0)**

One variable **increases**, the other **decreases**.

**Strong Negative (r = -0.8 to -1.0):**
```
Y
│ ●
│   ●
│     ●
│       ●
│         ●
└───────────── X

Points tightly clustered in downward direction
```

**Weak Negative (r = -0.3 to -0.7):**
```
Y
│●  ●
│  ●   ●
│    ●  ●
│      ●
└───────────── X

Points loosely scattered in downward trend
```

**Example:** Stress levels and sleep hours

---

#### **3. No Correlation (r ≈ 0)**

**Random scatter** with no discernible pattern.

```
Y
│ ●   ●   ●
│   ● ● ●
│ ●   ●   ●
│   ● ● ●
└───────────── X

Shotgun blast pattern - no relationship
```

**Example:** Shoe size and intelligence

---

### **2.1.3 Correlation Coefficient (r)**

**Pearson's Correlation Coefficient** quantifies the linear relationship.

**Computational Formula:**

```
        n(ΣXY) - (ΣX)(ΣY)
r = ─────────────────────────────────────────
    √{[nΣX² - (ΣX)²] × [nΣY² - (ΣY)²]}
```

Where:
- n = number of data pairs
- ΣXY = sum of products of paired X and Y values
- ΣX = sum of all X values
- ΣY = sum of all Y values
- ΣX² = sum of squared X values
- ΣY² = sum of squared Y values

---

### **2.1.4 Step-by-Step Calculation Examples**

#### **Example 1: Basic Correlation Calculation**

**Problem:** Calculate correlation coefficient for the following data:

| X | Y |
|---|---|
| 1 | 2 |
| 2 | 4 |
| 3 | 5 |
| 4 | 4 |
| 5 | 5 |

**Step-by-Step Solution:**

```
Step 1: Create calculation table

X    Y    X²    Y²    XY
─────────────────────────
1    2     1     4     2
2    4     4    16     8
3    5     9    25    15
4    4    16    16    16
5    5    25    25    25
─────────────────────────
ΣX=15 ΣY=20 ΣX²=55 ΣY²=86 ΣXY=66

n = 5

Step 2: Substitute into formula

        n(ΣXY) - (ΣX)(ΣY)
r = ─────────────────────────────────────────
    √{[nΣX² - (ΣX)²] × [nΣY² - (ΣY)²]}

        5(66) - (15)(20)
r = ────────────────────────────────────
    √{[5(55) - 15²] × [5(86) - 20²]}

        330 - 300
r = ─────────────────────────────
    √{[275 - 225] × [430 - 400]}

          30
r = ─────────────────────
    √{50 × 30}

        30           30
r = ───────── = ───────── = 0.7746
    √1500       38.73

Answer: r ≈ 0.77

Step 3: Interpretation

r = 0.77 indicates:
✓ Strong positive correlation
✓ As X increases, Y tends to increase
✓ Points cluster moderately tightly in upward direction
✓ About 60% (r² = 0.60) of variance in Y explained by X
```

---

#### **Example 2: Perfect Correlation**

**Problem:** Calculate r for:

| X | Y |
|---|---|
| 1 | 3 |
| 2 | 6 |
| 3 | 9 |
| 4 | 12 |

**Solution:**

```
X    Y    X²    Y²    XY
─────────────────────────
1    3     1     9     3
2    6     4    36    12
3    9     9    81    27
4   12    16   144    48
─────────────────────────
ΣX=10 ΣY=30 ΣX²=30 ΣY²=270 ΣXY=90

n = 4

        4(90) - (10)(30)
r = ────────────────────────────────────
    √{[4(30) - 10²] × [4(270) - 30²]}

        360 - 300
r = ──────────────────────────────
    √{[120 - 100] × [1080 - 900]}

         60
r = ─────────────────────
    √{20 × 180}

        60        60
r = ───────── = ──── = 1.0
    √3600       60

Answer: r = 1.0 (Perfect positive correlation)

Notice: Y = 3X (perfect linear relationship)
```

---

#### **Example 3: Negative Correlation**

**Problem:** Calculate r for study hours vs errors:

| Hours (X) | Errors (Y) |
|-----------|------------|
| 1 | 10 |
| 2 | 8 |
| 3 | 6 |
| 4 | 5 |
| 5 | 3 |

**Solution:**

```
X    Y    X²    Y²    XY
─────────────────────────
1   10     1   100    10
2    8     4    64    16
3    6     9    36    18
4    5    16    25    20
5    3    25     9    15
─────────────────────────
ΣX=15 ΣY=32 ΣX²=55 ΣY²=234 ΣXY=79

n = 5

        5(79) - (15)(32)
r = ────────────────────────────────────
    √{[5(55) - 15²] × [5(234) - 32²]}

        395 - 480
r = ───────────────────────────────
    √{[275 - 225] × [1170 - 1024]}

         -85
r = ──────────────────────
    √{50 × 146}

        -85         -85
r = ───────── = ───────── = -0.996
    √7300       85.44

Answer: r ≈ -0.996 (Very strong negative correlation)

Interpretation: More study hours strongly associated 
with fewer errors.
```

---

### **2.1.5 Interpreting r Values**

| r Value Range | Interpretation |
|---------------|----------------|
| **+0.8 to +1.0** | Strong positive correlation |
| **+0.3 to +0.7** | Moderate positive correlation |
| **0 to +0.2** | Weak or no positive correlation |
| **0** | No linear correlation |
| **-0.2 to 0** | Weak or no negative correlation |
| **-0.3 to -0.7** | Moderate negative correlation |
| **-0.8 to -1.0** | Strong negative correlation |

---

### **2.1.6 Coefficient of Determination (R²)**

**Definition:** The square of the correlation coefficient. Represents the **proportion of variance** in the dependent variable (Y) that is **explained** by the independent variable (X).

**Formula:**
```
R² = r²
```

**Interpretation:**
- If r = 0.9, then R² = 0.81
- This means **81%** of the variation in Y is explained by X
- The remaining **19%** is due to other unmeasured factors

**Example:**
```
r = 0.85
R² = 0.85² = 0.7225 = 72.25%

Interpretation: 72.25% of variance in Y is explained 
by X. The model has good explanatory power.
```

---

### **2.1.7 Pearson vs Spearman Correlation**

#### **Pearson Correlation (r)**

**Used for:**
- Continuous data
- Linear relationships
- Normally distributed data
- Interval or ratio scales

**Characteristics:**
- Measures strength of linear relationship
- Sensitive to outliers
- Parametric test

---

#### **Spearman Rank Correlation (ρ or rho)**

**Used for:**
- Ordinal (ranked) data
- Non-linear but monotonic relationships
- Non-normal distributions
- When outliers present

**Formula (Shortcut Method):**

```
          6Σd²
ρ = 1 - ───────────
        n(n² - 1)
```

Where:
- d = difference between ranks of X and Y
- n = number of pairs

---

### **Example: Spearman Rank Correlation**

**Problem:** Two judges rank 5 contestants:

| Contestant | Judge 1 (X) | Judge 2 (Y) |
|------------|-------------|-------------|
| A | 1 | 2 |
| B | 2 | 1 |
| C | 3 | 3 |
| D | 4 | 5 |
| E | 5 | 4 |

**Solution:**

```
Step 1: Create calculation table

Contestant   X   Y   d=X-Y   d²
────────────────────────────────
A            1   2    -1     1
B            2   1     1     1
C            3   3     0     0
D            4   5    -1     1
E            5   4     1     1
────────────────────────────────
                  Σd² = 4

n = 5

Step 2: Apply formula

          6Σd²
ρ = 1 - ───────────
        n(n² - 1)

          6(4)
ρ = 1 - ─────────────
        5(25 - 1)

          24         24
ρ = 1 - ───── = 1 - ──── = 1 - 0.2 = 0.8
        120        120

Answer: ρ = 0.8

Interpretation: Strong positive agreement between 
the two judges' rankings.
```

---

## 2.2 SCATTER PLOTS

### **2.2.1 Definition**

A **scatter plot** is a graphical representation that plots pairs of values (X, Y) as individual dots on a two-dimensional graph.

**Purpose:**
- Visualize the **direction** and **strength** of relationship
- Detect **outliers** visually
- Identify whether relationship is **linear or non-linear**
- Serve as starting point before computing r

---

### **2.2.2 Constructing a Scatter Plot**

**Steps:**

1. Place **independent variable (X)** on horizontal axis
2. Place **dependent variable (Y)** on vertical axis
3. Plot each (X, Y) pair as a dot
4. Observe the overall pattern

---

### **2.2.3 Reading Scatter Plots**

| Pattern | Interpretation |
|---------|----------------|
| Points tightly clustered in upward direction | Strong positive correlation |
| Points loosely scattered in upward direction | Weak positive correlation |
| Points randomly distributed (shotgun blast) | No correlation |
| Points tightly clustered in downward direction | Strong negative correlation |
| Points form a curve | Non-linear relationship |

---

### **2.2.4 Example: Scatter Plot Analysis**

**Data:** Study Hours vs Exam Scores

| Student | Hours (X) | Score (Y) |
|---------|-----------|-----------|
| 1 | 2 | 45 |
| 2 | 3 | 55 |
| 3 | 4 | 60 |
| 4 | 5 | 70 |
| 5 | 6 | 72 |
| 6 | 7 | 80 |
| 7 | 8 | 88 |
| 8 | 9 | 92 |

**ASCII Scatter Plot:**
```
Score (Y)
100 │                       ●
 95 │
 90 │                   ●
 85 │
 80 │               ●
 75 │
 70 │           ●
 65 │
 60 │       ●
 55 │   ●
 50 │
 45 │●
 40 │
    └─────────────────────────────
     0  1  2  3  4  5  6  7  8  9  10
           Study Hours (X)

Pattern: Strong positive correlation
As study hours increase, scores increase
Points cluster tightly along upward line
```

**Analysis:**
- Direction: Positive (upward trend)
- Strength: Strong (points tightly clustered)
- Form: Linear (straight line pattern)
- Outliers: None visible
- Expected r: Approximately +0.95 to +0.98

---

### **2.2.5 Relationship Between Correlation and Scatter Plots**

```
┌─────────────────────────────────────────┐
│         SCATTER PLOT                    │
│  ↓                                      │
│  Visual impression of relationship      │
│  ↓                                      │
│  Shows:                                 │
│  - Direction (positive/negative)        │
│  - Form (linear/non-linear)             │
│  - Strength (tight/loose clustering)    │
│  ↓                                      │
│  CALCULATE r                            │
│  ↓                                      │
│  Quantifies the relationship precisely  │
│  ↓                                      │
│  BOTH TOGETHER = Complete picture       │
└─────────────────────────────────────────┘
```

**Key Point:** 
- Scatter plot shows the **form**
- Correlation coefficient (r) quantifies the **strength**
- Always use BOTH for complete analysis

---

## 2.3 REGRESSION

### **2.3.1 Definition**

**Regression** is a method used to model the **functional relationship** between variables to **predict** a dependent variable (Y) based on an independent variable (X).

**Key Difference from Correlation:**
- **Correlation**: Measures strength (symmetric, no predictions)
- **Regression**: Models predictive functional relationships (asymmetric, X predicts Y)

---

### **2.3.2 Regression Line (Line of Best Fit)**

**Definition:** A straight line that **best summarizes** the relationship between two variables.

**Equation:**
```
Ŷ = a + bX
```

Where:
- **Ŷ** (Y-hat) = Predicted value of Y
- **a** = Y-intercept (predicted Y when X = 0)
- **b** = Slope (change in Y for every one-unit increase in X)
- **X** = Independent variable

---

### **2.3.3 Least Squares Principle**

**Definition:** The mathematical method used to determine the best-fitting regression line.

**Concept:**
Finds the line that **minimizes the sum of squared vertical distances (residuals)** between each observed data point and the corresponding point on the regression line.

**Mathematical Expression:**
```
Minimize Σ(Y - Ŷ)²
```

Where:
- Y = Actual observed value
- Ŷ = Predicted value from regression line
- (Y - Ŷ) = Residual (error)

**Why Squared?**
- Eliminates negative signs
- Penalizes larger errors more heavily
- Mathematically tractable

---

### **2.3.4 Calculating Regression Coefficients**

**Formula for Slope (b):**

```
        n(ΣXY) - (ΣX)(ΣY)
b = ─────────────────────────
        nΣX² - (ΣX)²
```

**Formula for Intercept (a):**

```
a = Ȳ - bX̄
```

Where:
- X̄ = Mean of X
- Ȳ = Mean of Y

---

### **2.3.5 Step-by-Step Regression Example**

#### **Example 1: Complete Regression Analysis**

**Problem:** Develop regression equation to predict exam scores (Y) from study hours (X).

**Data:**

| Student | Hours (X) | Score (Y) |
|---------|-----------|-----------|
| 1 | 1 | 50 |
| 2 | 2 | 55 |
| 3 | 3 | 65 |
| 4 | 4 | 70 |
| 5 | 5 | 80 |

**Step-by-Step Solution:**

```
Step 1: Create calculation table

X    Y    X²    Y²    XY
─────────────────────────
1   50     1  2500    50
2   55     4  3025   110
3   65     9  4225   195
4   70    16  4900   280
5   80    25  6400   400
─────────────────────────
ΣX=15 ΣY=320 ΣX²=55 ΣY²=21050 ΣXY=1035

n = 5

Step 2: Calculate means

X̄ = ΣX/n = 15/5 = 3
Ȳ = ΣY/n = 320/5 = 64

Step 3: Calculate slope (b)

        n(ΣXY) - (ΣX)(ΣY)
b = ─────────────────────────
        nΣX² - (ΣX)²

        5(1035) - (15)(320)
b = ───────────────────────────
        5(55) - 15²

        5175 - 4800
b = ─────────────────
        275 - 225

         375
b = ───────── = 7.5
         50

Interpretation: Each additional study hour 
increases score by 7.5 points.

Step 4: Calculate intercept (a)

a = Ȳ - bX̄
a = 64 - 7.5(3)
a = 64 - 22.5
a = 41.5

Interpretation: With zero study hours, 
predicted score is 41.5.

Step 5: Write regression equation

Ŷ = 41.5 + 7.5X

Step 6: Make predictions

Question: Predict score for 6 hours of study

Ŷ = 41.5 + 7.5(6)
Ŷ = 41.5 + 45
Ŷ = 86.5

Answer: Predicted score = 86.5

Question: Predict score for 3.5 hours

Ŷ = 41.5 + 7.5(3.5)
Ŷ = 41.5 + 26.25
Ŷ = 67.75

Answer: Predicted score = 67.75

Step 7: Visualize regression line

Score (Y)
 90 │                       ● (6, 86.5)
    │                     ╱
 85 │                   ╱
    │                 ╱
 80 │               ● (5, 80)
    │             ╱
 75 │           ╱
    │         ╱
 70 │       ● (4, 70)
    │     ╱
 65 │   ● (3, 65)
    │ ╱
 60 │╱
    │
 55 │  ● (2, 55)
    │
 50 │● (1, 50)
    │
 45 │
    └─────────────────────────
     0  1  2  3  4  5  6  7
         Study Hours (X)

Regression line: Ŷ = 41.5 + 7.5X
```

---

#### **Example 2: Real-World Application**

**Scenario:** A retail company wants to predict sales based on advertising spend.

**Data (Monthly):**

| Month | Ad Spend (₹K) | Sales (₹L) |
|-------|---------------|------------|
| Jan | 10 | 120 |
| Feb | 15 | 150 |
| Mar | 20 | 180 |
| Apr | 25 | 210 |
| May | 30 | 250 |

**Solution:**

```
Step 1: Calculation table

X     Y     X²     Y²      XY
──────────────────────────────
10   120   100   14400    1200
15   150   225   22500    2250
20   180   400   32400    3600
25   210   625   44100    5250
30   250   900   62500    7500
──────────────────────────────
ΣX=100 ΣY=910 ΣX²=2250 ΣY²=175900 ΣXY=19800

n = 5

Step 2: Means

X̄ = 100/5 = 20
Ȳ = 910/5 = 182

Step 3: Slope

        5(19800) - (100)(910)
b = ───────────────────────────
        5(2250) - 100²

        99000 - 91000
b = ───────────────────
        11250 - 10000

         8000
b = ───────── = 6.4
         1250

Interpretation: Each ₹1K increase in ad spend 
generates ₹6.4L additional sales.

Step 4: Intercept

a = 182 - 6.4(20)
a = 182 - 128
a = 54

Interpretation: Base sales without advertising = ₹54L

Step 5: Regression equation

Ŷ = 54 + 6.4X

Step 6: Business application

If company spends ₹35K on ads:

Ŷ = 54 + 6.4(35)
Ŷ = 54 + 224
Ŷ = ₹278L

Predicted sales = ₹278 lakhs

ROI Analysis:
Ad spend: ₹35K = ₹0.35L
Sales: ₹278L
Return ratio: 278/0.35 = 794:1

Excellent ROI!
```

---

### **2.3.6 Types of Regression**

#### **1. Simple Linear Regression**
- One independent variable (X) predicts one dependent variable (Y)
- Equation: Ŷ = a + bX
- Example: Predicting house price from size only

#### **2. Multiple Linear Regression**
- Two or more independent variables predict one dependent variable
- Equation: Ŷ = β₀ + β₁X₁ + β₂X₂ + ... + βₙXₙ
- Example: Predicting house price from size, bedrooms, age, location

#### **3. Linear Regression**
- Assumes relationship forms a straight line
- Most common type

#### **4. Non-linear Regression**
- Relationship forms a curve
- Examples: Quadratic, exponential, logarithmic

---

## 2.4 MULTIPLE REGRESSION

### **2.4.1 Definition**

**Multiple Linear Regression** is an extension of simple linear regression that models the relationship between **one dependent variable (Y)** and **two or more independent variables (X₁, X₂, ... Xₙ)**.

It fits a **hyperplane in n-dimensional space** that best predicts Y from the combined set of predictors.

---

### **2.4.2 The Regression Equation**

```
Ŷ = β₀ + β₁X₁ + β₂X₂ + β₃X₃ + ... + βₙXₙ + ε
```

**Term Definitions:**

| Term | Name | Meaning |
|------|------|---------|
| **Ŷ** | Predicted value | Model's best estimate of Y |
| **β₀** | Intercept | Value of Y when ALL predictors = 0 |
| **β₁, β₂...** | Coefficients/Slopes | Change in Y for 1-unit increase in Xᵢ, holding all others constant |
| **Xᵢ** | Predictors | Independent/input variables |
| **ε** | Error term | Unexplained residual, assumed normally distributed |

---

### **2.4.3 How Regression Coefficients Indicate Relationships**

**Key Concept:** Each coefficient βᵢ gives the **isolated, partial effect** of one predictor on Y while **statistically controlling for all other variables**.

This is the **core power of multiple regression**.

---

### **2.4.4 Worked Example: Predicting House Price**

**Scenario:** Predict house price (in ₹ lakhs) using:
- Size (m²)
- Age (years)
- Number of bedrooms

**Regression Equation:**
```
Price = 10 + 0.8(Size) - 0.5(Age) + 3.2(Bedrooms)
```

**Interpretation of Coefficients:**

```
β₀ = 10
→ Base price of ₹10L when all inputs are theoretically zero

β₁ = +0.8 (Size)
→ Each additional square meter adds ₹0.8L to price
→ Holding age and bedrooms constant
→ Positive relationship

β₂ = -0.5 (Age)
→ Each additional year of age reduces price by ₹0.5L
→ Holding other variables constant
→ Negative relationship

β₃ = +3.2 (Bedrooms)
→ Each additional bedroom adds ₹3.2L
→ Holding other variables constant
→ Positive relationship
```

**Prediction Example:**

**Problem:** Predict price for a 120 m² house, 5 years old, with 3 bedrooms.

**Solution:**

```
Price = 10 + 0.8(Size) - 0.5(Age) + 3.2(Bedrooms)

Price = 10 + 0.8(120) - 0.5(5) + 3.2(3)

Price = 10 + 96 - 2.5 + 9.6

Price = ₹113.1 Lakhs

Answer: Predicted house price = ₹113.1 lakhs
```

**Key Insight:**
Each coefficient tells you the **isolated effect** of one variable while all others stay constant — this is the power of multiple regression.

---

### **2.4.5 Key Assumptions of Multiple Regression**

For valid results, these **six key assumptions** must hold:

#### **1. Linearity**
- Relationship between each predictor and Y must be linear
- **Check:** Scatter plots, residual vs. fitted plots
- **Violation:** Transform variables or use non-linear models

#### **2. No Multicollinearity**
- Predictors must NOT be highly correlated with each other
- **Check:** Variance Inflation Factor (VIF)
  - VIF < 10: Acceptable
  - VIF > 10: Problematic multicollinearity
- **Solution:** Remove or combine highly correlated predictors

#### **3. Homoscedasticity**
- Residuals have **constant variance** across all predicted values
- **Check:** Residuals vs. Fitted plot
  - Should show random scatter (no funnel shape)
- **Violation:** Use robust regression or transform variables

#### **4. Independence of Errors**
- Residuals are independent of each other
- **Critical for:** Time-series data
- **Check:** Durbin-Watson test
- **Violation:** Use time-series models (ARIMA, etc.)

#### **5. Normality of Residuals**
- Errors follow a normal distribution
- **Check:** 
  - Q-Q plot (points should follow diagonal line)
  - Shapiro-Wilk test
- **Violation:** Transform dependent variable

#### **6. No Endogeneity**
- Predictors must NOT be correlated with the error term
- **Violated by:** Omitted variable bias
- **Solution:** Include all relevant variables

---

### **2.4.6 Advantages Over Simple Regression**

1. **Explains more variance in Y**
   - Incorporates multiple factors simultaneously
   - Higher R² value

2. **Controls for confounding variables**
   - Gives purer estimates of each predictor's effect
   - Isolates individual contributions

3. **Enables more accurate predictions**
   - Realistic modeling of complex scenarios
   - Better generalization

---

### **2.4.7 Comprehensive Retail Sales Example**

**Scenario:** A retail company wants to predict monthly sales (Y) using:
- TV advertising spend (X₁)
- Social media advertising spend (X₂)
- Store size in sq ft (X₃)
- Location type: Urban=1, Rural=0 (X₄)

**Step 1: Define the Model**

```
Ŷ(Sales) = β₀ + β₁(TV Ads) + β₂(Social Media) + β₃(Store Size) + β₄(Location)
```

**Step 2: Hypothetical Data Collection**

| Month | TV Ads (₹K) | Social Media (₹K) | Store Size (sq ft) | Location (1=Urban) | Sales (₹L) |
|-------|-------------|-------------------|-------------------|-------------------|------------|
| Jan | 50 | 20 | 1500 | 1 | 120 |
| Feb | 30 | 10 | 1200 | 0 | 75 |
| Mar | 70 | 35 | 2000 | 1 | 180 |
| Apr | 45 | 25 | 1800 | 1 | 145 |
| May | 20 | 8 | 1000 | 0 | 60 |

**Step 3: Build the Model (Python Code)**

```python
import pandas as pd
from sklearn.linear_model import LinearRegression

data = {
    'TV_Ads': [50, 30, 70, 45, 20],
    'Social_Media': [20, 10, 35, 25, 8],
    'Store_Size': [1500, 1200, 2000, 1800, 1000],
    'Location': [1, 0, 1, 1, 0],
    'Sales': [120, 75, 180, 145, 60]
}

df = pd.DataFrame(data)

X = df[['TV_Ads', 'Social_Media', 'Store_Size', 'Location']]
Y = df['Sales']

model = LinearRegression()
model.fit(X, Y)

print("Intercept (β₀):", model.intercept_)
print("Coefficients:", dict(zip(X.columns, model.coef_)))
print("R² Score:", model.score(X, Y))
```

**Step 4: Hypothetical Output and Interpretation**

Suppose the model produces:

```
Sales = 5 + 1.8(TV Ads) + 0.9(Social Media) + 0.04(Store Size) + 12(Location)
```

**Coefficient Interpretation:**

| Coefficient | Value | Interpretation |
|-------------|-------|----------------|
| **β₀** | 5 | Base sales of ₹5L when all inputs are zero |
| **β₁** | 1.8 | Each ₹1K increase in TV ads → +₹1.8L in sales |
| **β₂** | 0.9 | Each ₹1K in social media → +₹0.9L in sales |
| **β₃** | 0.04 | Each additional sq ft → +₹0.04L in sales |
| **β₄** | 12 | Urban stores sell ₹12L more than rural stores |

**Step 5: Predict Future Sales**

**Problem:** Predict sales for a new store with:
- TV Ads = ₹60K
- Social Media = ₹30K
- Store Size = 1700 sq ft
- Urban location

**Solution:**

```
Sales = 5 + 1.8(TV Ads) + 0.9(Social Media) + 0.04(Store Size) + 12(Location)

Sales = 5 + 1.8(60) + 0.9(30) + 0.04(1700) + 12(1)

Sales = 5 + 108 + 27 + 68 + 12

Sales = ₹220 Lakhs

Answer: Predicted sales = ₹220 lakhs
```

**Step 6: Analysis of Key Drivers**

```
✓ TV Ads (β = 1.8) is the strongest advertising driver
  → Highest return per rupee spent

✓ Location (β = 12) has significant fixed advantage
  → Urban stores perform better regardless of other factors

✓ Social Media (β = 0.9) contributes positively
  → But less effective than TV

✓ Store Size (β = 0.04) has modest effect
  → Larger stores generate slightly more sales

MANAGEMENT RECOMMENDATIONS:
1. Prioritize TV advertising budget
2. Expand in urban locations
3. Social media is secondary priority
4. Store size optimization less critical

If R² close to 1:
→ Model explains most variance in sales
→ High predictive accuracy
```

---

## 2.5 STANDARD ERROR OF ESTIMATE

### **2.5.1 Definition**

The **Standard Error of Estimate (SEE)** measures **prediction accuracy** by calculating the **average distance** between the observed values and the predicted values on the regression line.

**Formula:**
```
         _____________
        / Σ(Y - Ŷ)²
SEE = √ ────────────
           n - 2
```

Where:
- Y = Actual observed value
- Ŷ = Predicted value
- (Y - Ŷ) = Residual
- n = Sample size

---

### **2.5.2 Interpretation**

- **Smaller SEE** → Regression line fits data better → Predictions more accurate
- **Larger SEE** → More scatter around line → Predictions less accurate
- Expressed in **same units** as the Y variable
- **Higher R²** (more variance explained) → **Lower SEE**

---

### **2.5.3 Example Calculation**

**Problem:** Calculate SEE for the following regression:

| X | Y | Ŷ (Predicted) |
|---|---|---------------|
| 1 | 50 | 49 |
| 2 | 55 | 56.5 |
| 3 | 65 | 64 |
| 4 | 70 | 71.5 |
| 5 | 80 | 79 |

**Solution:**

```
Step 1: Calculate residuals (Y - Ŷ)

X   Y   Ŷ    Y - Ŷ    (Y - Ŷ)²
───────────────────────────────
1  50  49     1.0       1.00
2  55  56.5  -1.5       2.25
3  65  64     1.0       1.00
4  70  71.5  -1.5       2.25
5  80  79     1.0       1.00
───────────────────────────────
                Σ(Y-Ŷ)² = 7.50

Step 2: Apply formula

         _____________
        / Σ(Y - Ŷ)²
SEE = √ ─────────────
           n - 2

         _______
        / 7.50
SEE = √ ───────
          5-2

        _____
SEE = √2.5 = 1.58

Interpretation:
On average, actual scores deviate from predicted 
scores by about 1.58 points.

This is relatively small compared to the score 
range (50-80), indicating good model fit.
```

---

## 2.6 REGRESSION TO THE MEAN

### **2.6.1 Definition**

**Regression to the Mean (RTM)** is the statistical tendency for **extreme measurements** to move **closer to the average** on a subsequent measurement, purely due to **random variation** — NOT because of any real change or intervention.

---

### **2.6.2 Why It Happens**

**Concept:**
Any measurement contains:
- **True score** (actual ability/characteristic)
- **Random error** (luck, measurement error, temporary factors)

**Extreme scores** have unusually large random errors that do not persist.

```
┌─────────────────────────────────────────┐
│         MEASUREMENT COMPOSITION         │
├─────────────────────────────────────────┤
│                                         │
│  Observed Score = True Score + Error   │
│                                         │
│  Example:                              │
│  Exam Score = Actual Knowledge + Luck  │
│                                         │
│  Extreme HIGH score:                   │
│  → High knowledge + Unusually good luck│
│                                         │
│  Next attempt:                         │
│  → High knowledge + Normal luck        │
│  → Score moves closer to true ability  │
└─────────────────────────────────────────┘
```

---

### **2.6.3 Mathematical Formula**

```
E[Y₂ | Y₁] = μ + ρ × (Y₁ - μ)
```

Where:
- **E[Y₂ | Y₁]** = Expected 2nd measurement given the first
- **μ (mu)** = Population mean (true average)
- **ρ (rho)** = Correlation between measurements (0 to 1)
- **Y₁ - μ** = Deviation from mean on first measurement

**Key Insight:**
- The further Y₁ is from the mean AND
- The weaker the correlation (ρ)
- → The stronger the regression to the mean effect

**Special Cases:**
- If ρ = 1 (perfect reliability) → No RTM
- If ρ = 0 (no reliability) → Full RTM to mean
- In practice: 0 < ρ < 1

---

### **2.6.4 Visualizing Regression to the Mean**

```
Attempt 1 Scores          Attempt 2 Scores (regressed)

     ● Extreme HIGH  →         ● Closer to mean
                               │
                               │
         ●                     ●
                               │
         │                     │
    ─────┼───── Mean ──────────┼─────
         │                     │
                               │
         ●                     ●
                               │
                               │
     ● Extreme LOW   →         ● Closer to mean

What's Happening:
✓ Extreme HIGH scorers tend to score LOWER next time
✓ Extreme LOW scorers tend to score HIGHER next time
✓ Average scorers stay near average
✓ Purely due to random variation - NOT improvement/decline!
```

---

### **2.6.5 Real-World Examples**

#### **Example 1: Sports — Sophomore Slump**

**Scenario:**
A cricketer has a record-breaking debut season with exceptional performance (e.g., 1000 runs).

**Next Season:**
Performance drops to 650 runs (closer to career average).

**Common Mistake:**
Analysts blame "extra pressure" or "opponents figuring him out."

**RTM Explanation:**
- Debut season: True ability + Unusually good luck
- Second season: True ability + Normal luck
- The drop was **statistical, not causal**

**⚠️ Warning:** Don't attribute RTM to interventions without control groups!

---

#### **Example 2: Medicine — Patient Blood Pressure**

**Scenario:**
Patients with dangerously high BP (e.g., 180/120) are enrolled in a study.

**After 2 weeks:**
BP drops to 150/95 — even in the placebo group!

**Common Mistake:**
Doctors credit the treatment.

**RTM Explanation:**
- Initial extreme reading: True BP + Measurement error/stress
- Follow-up reading: True BP + Normal variation
- Patients with extreme readings naturally have more moderate readings on follow-up

**Solution:** Use control groups to isolate true treatment effect!

---

#### **Example 3: Education — Tutoring Studies**

**Scenario:**
Students who scored in the bottom 10% are enrolled in a tutoring program.

**After tutoring:**
Scores improve significantly.

**Question:** Is it the tutoring?

**RTM Concern:**
- Initial low score: True ability + Bad luck/anxiety
- Next test: True ability + Normal conditions
- Part of the gain is RTM

**Solution:** Need a control group (similar low-scoring students without tutoring) to isolate the true tutoring effect.

---

#### **Example 4: Finance — Fund Manager Performance**

**Scenario:**
Top-performing mutual funds in Year 1 often perform near-average in Year 2.

**Common Mistake:**
Investors assume a strategy shift is needed.

**RTM Explanation:**
- Year 1: Skill + Exceptional luck
- Year 2: Skill + Normal luck
- The "curse of the top performer"

**⚠️ Warning:** This is one of the most common statistical errors in finance!

---

### **2.6.6 How to Avoid RTM Fallacy**

**Common Pitfall:**
Attributing RTM to an intervention

**Solution:**
✅ Always include a **control group**
✅ Use **randomized controlled trials**
✅ Natural regression to the mean can look like a treatment effect in pre-post studies

**Example of Proper Design:**

```
┌─────────────────────────────────────────┐
│         PROPER EXPERIMENTAL DESIGN      │
├─────────────────────────────────────────┤
│                                         │
│  Treatment Group:                      │
│  Low scorers → Tutoring → Test again  │
│                                         │
│  Control Group:                        │
│  Low scorers → No tutoring → Test again│
│                                         │
│  Compare:                              │
│  (Treatment improvement) -             │
│  (Control improvement)                 │
│  = TRUE tutoring effect                │
│                                         │
│  The control group accounts for RTM!   │
└─────────────────────────────────────────┘
```

---

### **2.6.7 Francis Galton's Discovery (1886)**

**Historical Context:**
First described while studying heights of parents and children.

**Finding:**
- Tall parents have tall children
- BUT children are closer to average height than the parents
- Short parents have short children
- BUT children are closer to average height than the parents

**Conclusion:**
Extreme values "regress" toward the mean across generations.

---

## 2.7 BASICS OF NUMPY ARRAYS

### **2.7.1 Introduction to NumPy**

**NumPy (Numerical Python)** is the core library for scientific computing in Python.

**Key Features:**
- Provides powerful **N-dimensional array object (ndarray)**
- Designed for basic and advanced array operations
- Much faster than standard Python lists
- Foundation for all data science libraries

---

### **2.7.2 Why NumPy Arrays Over Python Lists?**

| Feature | NumPy Array | Python List |
|---------|-------------|-------------|
| **Data Type** | Homogeneous (single type) | Heterogeneous (mixed types) |
| **Memory** | Compact, efficient | Extra overhead per element |
| **Speed** | Orders of magnitude faster | Slower (Python loops) |
| **Operations** | Vectorized (element-wise) | Requires explicit loops |
| **Functionality** | Built-in math functions | Limited |

---

### **2.7.3 Creating NumPy Arrays**

```python
import numpy as np

# From a list
arr1 = np.array([1, 2, 3, 4, 5])

# 2D array (matrix)
arr2 = np.array([[1, 2, 3], 
                 [4, 5, 6]])

# Array of zeros
zeros = np.zeros((3, 4))  # 3 rows, 4 columns

# Array of ones
ones = np.ones((2, 3))

# Array with range
range_arr = np.arange(0, 10, 2)  # [0, 2, 4, 6, 8]

# Evenly spaced values
linspace_arr = np.linspace(0, 1, 5)  # [0, 0.25, 0.5, 0.75, 1]

# Random numbers
random_arr = np.random.randn(3, 3)  # 3x3 normal distribution
```

---

### **2.7.4 NumPy Array Attributes**

**Example:**

```python
import numpy as np

np.random.seed(0)  # Initialize random number generator

x1 = np.random.randint(10, size=5)       # 1D array
x2 = np.random.randint(10, size=(2, 4))  # 2D array
x3 = np.random.randint(10, size=(3, 4, 5))  # 3D array

print("x3 ndim:", x3.ndim)      # → 3 (Number of dimensions)
print("x3 shape:", x3.shape)    # → (3, 4, 5) (Size of each dimension)
print("x3 size:", x3.size)      # → 60 [3×4×5] (Total elements)
print("x3 dtype:", x3.dtype)    # → int64 (Data type)
print("itemsize:", x3.itemsize) # → 8 bytes (Size per element)
print("nbytes:", x3.nbytes)     # → 480 bytes [60×8] (Total size)
```

**ASCII Visualization of Dimensions:**

```
1D Array          2D Array (Matrix)      3D Array (Depth)
+---+            +---+---+---+---+       /---+---+---+---+
| 0 |            |0,0|0,1|0,2|0,3|      / / / / /
+---+            +---+---+---+---+     / / / / /
| 1 |            |1,0|1,1|1,2|1,3|    / / / / /
+---+            +---+---+---+---+   / / / / /
| 2 |                              / / / / /
+---+                             / / / / /
                                 / / / / /
                                +---+---+---+---+
```

---

### **2.7.5 Array Indexing - Accessing Single Elements**

```python
import numpy as np

# 1D indexing
x1 = np.array([5, 0, 3, 3, 7, 9])

x1[0]   # → 5 (first element)
x1[-1]  # → 9 (last element, negative index)
x1[-2]  # → 7 (second from last)

# 2D indexing (comma-separated tuple)
x2 = np.array([[3, 5, 2, 4],
               [7, 6, 8, 8],
               [1, 6, 7, 7]])

x2[0, 0]   # → 3 (row 0, col 0)
x2[2, 0]   # → 1 (row 2, col 0)
x2[2, -1]  # → 7 (row 2, last col)

# 3D indexing
arr = np.array([[[1, 2, 3], [4, 5, 6]],
                [[7, 8, 9], [10, 11, 12]]])

print(arr[0, 1, 2])  # → 6
# Accessing: 3rd element of 2nd array of 1st array
```

---

### **2.7.6 Array Slicing - Accessing Subarrays**

**Syntax:** `x[start:stop:step]`

**Examples:**

```python
import numpy as np

x = np.arange(10)  # [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

x[:5]      # → [0, 1, 2, 3, 4] (first 5)
x[5:]      # → [5, 6, 7, 8, 9] (from index 5)
x[4:7]     # → [4, 5, 6] (middle subarray)
x[::2]     # → [0, 2, 4, 6, 8] (every other element)
x[1::2]    # → [1, 3, 5, 7, 9] (every other, starting from 1)
x[::-1]    # → [9, 8, 7, 6, 5, 4, 3, 2, 1, 0] (reversed)

# Multi-dimensional slicing
x2 = np.array([[20, 5, 2, 4],
               [7, 6, 8, 8],
               [1, 6, 7, 7]])

x2[:2, :3]    # First 2 rows, first 3 cols
x2[:3, ::2]   # All rows, every other column
x2[:, 0]      # All rows, col 0 → [20, 7, 1]
x2[0, :]      # Row 0, all cols → [20, 5, 2, 4]
```

---

### **2.7.7 Reshaping Arrays**

```python
import numpy as np

# Reshape 1D to 2D
grid = np.arange(1, 10).reshape((3, 3))
# → [[1, 2, 3],
#    [4, 5, 6],
#    [7, 8, 9]]

# Row vector
x = np.array([1, 2, 3])
x.reshape((1, 3))  # → [[1, 2, 3]]

# Column vector
x.reshape((3, 1))  # → [[1],
                   #    [2],
                   #    [3]]
```

---

### **2.7.8 Concatenation and Splitting**

**Concatenation:**

```python
import numpy as np

x = np.array([1, 2, 3])
y = np.array([3, 2, 1])

# 1D concatenation
np.concatenate([x, y])  # → [1, 2, 3, 3, 2, 1]

# 2D concatenation
grid2 = np.array([[1, 2, 3],
                  [4, 5, 6]])

np.concatenate([grid2, grid2])        # Vertical (axis=0)
np.concatenate([grid2, grid2], axis=1)  # Horizontal (axis=1)

# Stacking
np.vstack([x, grid2])  # Vertical stack (cols must match)
np.hstack([grid2, grid2])  # Horizontal stack (rows must match)
```

**Splitting:**

```python
arr = np.arange(8)  # [0, 1, 2, 3, 4, 5, 6, 7]

x1, x2, x3 = np.split(arr, [3, 5])
# x1 = [0, 1, 2]
# x2 = [3, 4]
# x3 = [5, 6, 7]

# 2D splitting
np.vsplit(grid2, 2)  # Vertical split
np.hsplit(grid2, 3)  # Horizontal split
```

---

## 2.8 NUMPY AGGREGATIONS & UNIVERSAL FUNCTIONS

### **2.8.1 Aggregation Functions**

Aggregations summarize an entire array (or along an axis) into a single (or fewer) values.

**Examples:**

```python
import numpy as np

L = np.random.random(100)  # 100 random numbers

# Basic aggregations
np.sum(L)     # Sum of all elements
np.min(L)     # Minimum value
np.max(L)     # Maximum value
np.mean(L)    # Mean (average)
np.std(L)     # Standard deviation
np.var(L)     # Variance
np.median(L)  # Median

# Product
np.prod(L)    # Product of all elements

# Index of min/max
np.argmin(L)  # Index of minimum value
np.argmax(L)  # Index of maximum value

# Percentile
np.percentile(L, 25)  # 25th percentile

# Boolean aggregations
np.any(L > 0.5)  # True if ANY element > 0.5
np.all(L > 0)    # True if ALL elements > 0

# Multi-dimensional aggregations
M = np.random.random((3, 4))

M.min(axis=0)  # Min of each column
M.max(axis=1)  # Max of each row
M.mean(axis=0) # Mean of each column
```

**Key Difference:**
- **Aggregations** reduce dimensions
- **Computations (UFuncs)** maintain same shape

---

### **2.8.2 Universal Functions (UFuncs)**

A **universal function (ufunc)** operates on ndarrays in an **element-by-element** fashion.

**Arithmetic UFuncs:**

```python
import numpy as np

x = np.array([1, 2, 3, 4])

x + 5          # np.add(x, 5)       → [6, 7, 8, 9]
x - 2          # np.subtract(x, 2)  → [-1, 0, 1, 2]
x * 3          # np.multiply(x, 3)  → [3, 6, 9, 12]
x / 2          # np.divide(x, 2)    → [0.5, 1.0, 1.5, 2.0]
x // 2         # np.floor_divide()  → [0, 1, 1, 2]
x ** 2         # np.power(x, 2)     → [1, 4, 9, 16]
x % 2          # np.mod(x, 2)       → [1, 0, 1, 0]
-x             # np.negative(x)     → [-1, -2, -3, -4]
```

---

### **2.8.3 Trigonometric Functions**

**Important:** Angles must be in **RADIANS**!

```python
import numpy as np

# Convert degrees to radians
Arr = np.array([0, 30, 60, 90])
arr = Arr * np.pi / 180  # Convert to radians

# Basic trig
np.sin(arr)   # → [0.0, 0.5, 0.866, 1.0]
np.cos(arr)   # → [1.0, 0.866, 0.5, 0.0]
np.tan(arr)   # → [0.0, 0.577, 1.732, ∞]

# Inverse trig
np.arcsin(arr)
np.arccos(arr)
np.arctan(arr)

# Hyperbolic
np.sinh(arr)
np.cosh(arr)
np.tanh(arr)

# Conversion
np.deg2rad(Arr)  # Degrees → Radians
np.rad2deg(arr)  # Radians → Degrees
```

---

### **2.8.4 Exponents and Logarithms**

```python
import numpy as np

x = np.array([1, 2, 4, 10])

# Exponents
np.exp(x)        # e^x
np.exp2(x)       # 2^x
np.power(3, x)   # 3^x

# Logarithms
np.log(x)    # Natural log (ln)
np.log2(x)   # Base-2 log
np.log10(x)  # Base-10 log
```

---

### **2.8.5 Other Mathematical Functions**

```python
import numpy as np

# Rounding
arr = np.array([3.14159, 2.71828])
np.around(arr, 3)  # Round to 3 decimals → [3.142, 2.718]
np.floor(arr)      # Round down → [3., 2.]
np.ceil(arr)       # Round up → [4., 3.]

# Absolute values
np.abs(np.array([-2, -1, 0, 1, 2]))  # → [2, 1, 0, 1, 2]
```

---

## 2.9 COMPARISONS, MASKS & BOOLEAN LOGIC

### **2.9.1 Comparison Operators**

Comparison operators return **boolean arrays** (True/False per element).

```python
import numpy as np

a = np.reshape(np.arange(16), (4, 4))
# → [[ 0,  1,  2,  3],
#    [ 4,  5,  6,  7],
#    [ 8,  9, 10, 11],
#    [12, 13, 14, 15]]

# Comparison returns boolean array
large_values = (a > 10)
# → [[False, False, False, False],
#    [False, False, False, False],
#    [False, False, False,  True],
#    [ True,  True,  True,  True]]

# Boolean operators
np.equal(a, 5)        # == element-wise
np.not_equal(a, 5)    # !=
np.less(a, 5)         # <
np.less_equal(a, 5)   # <=
np.greater(a, 5)      # >
np.greater_equal(a, 5) # >=
```

---

### **2.9.2 Logical Operators on Boolean Arrays**

```python
import numpy as np

a = np.arange(16)

# NOT (negation)
~(a > 10)

# AND
(a > 5) & (a < 12)

# OR
(a < 5) | (a > 12)

# XOR
(a < 5) ^ (a > 12)
```

---

### **2.9.3 Boolean Masking**

**Masking** means extracting, modifying, counting, or manipulating values based on some criterion.

**Examples:**

```python
import numpy as np

data = np.array([15, 42, 8, 73, 29, 55, 3, 90])

# Create a boolean mask
mask = data > 30
# mask → [False, True, False, True, False, True, False, True]

# Apply the mask to filter values
filtered = data[mask]
# filtered → [42, 73, 55, 90]

# Combined condition (values between 20 and 60)
filtered2 = data[(data > 20) & (data < 60)]
# filtered2 → [42, 29, 55]

# Count elements meeting condition
count = np.sum(data > 30)  # → 4

# Extract all even elements
even_values = (a % 2 == 0)
print(a[even_values])  # All even elements
```

**Use Case:** Extremely useful in data analysis for selecting rows meeting specific criteria without writing explicit loops.

---

## 2.10 FANCY INDEXING & STRUCTURED ARRAYS

### **2.10.1 Fancy Indexing**

**Definition:** Indexing an array with another NumPy array, a Python list, or a sequence of integers whose values select elements.

**Examples:**

```python
import numpy as np

A = np.linspace(0, 1, 11)
# → [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# Using another array to index
print(A[np.array([0, 2, 4])])  # → [0.0, 0.2, 0.4]

# Using a list to index
print(A[[0, 2, 4]])  # → [0.0, 0.2, 0.4]

# More complex example
arr = np.array([10, 20, 30, 40, 50, 60, 70, 80])

# Select elements at indices 0, 3, and 6
print(arr[[0, 3, 6]])  # → [10, 40, 70]

# Using a NumPy index array
indices = np.array([1, 4, 7])
print(arr[indices])  # → [20, 50, 80]

# 2D fancy indexing
B = np.array([[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]])

print(B[[0, 2], [1, 2]])  # → [2, 9]
# (row 0 col 1, row 2 col 2)
```

**Useful in Data Analysis For:**
- Selecting specific rows/records of interest
- Reordering or shuffling arrays (e.g., creating train/test splits)
- Extracting non-contiguous elements that cannot be captured by simple slicing

---

### **2.10.2 Structured Arrays**

**Definition:** Unlike regular NumPy arrays (which are homogeneous), structured arrays allow **different data types in each column** using a dtype definition.

**Method 1: Dictionary**

```python
import numpy as np

dtype = np.dtype({
    'names': ('name', 'age', 'weight'),
    'formats': ('U10', 'i4', 'f8')
})

data = np.array([('Alice', 25, 55.0),
                 ('Bob', 30, 75.5)], dtype=dtype)

print(data['name'])   # → ['Alice' 'Bob']
print(data['age'])    # → [25 30]
print(data['weight']) # → [55.0 75.5]
```

**Method 2: List of Tuples**

```python
dtype2 = np.dtype([('name', 'U10'),
                   ('age', 'i4'),
                   ('weight', 'f8')])

data2 = np.array([('Alice', 25, 55.0),
                  ('Bob', 30, 75.5)], dtype=dtype2)
```

**NumPy Type Codes:**

| Code | Type |
|------|------|
| **U** | Unicode string |
| **i** | Signed integer |
| **f** | Float |
| **b** | Boolean |
| **u** | Unsigned integer |
| **c** | Complex |
| **m** | Timedelta |
| **M** | Datetime |

---

## 2.11 DATA MANIPULATION WITH PANDAS

### **2.11.1 Introduction to Pandas**

**Pandas** is a high-level data manipulation tool built on the NumPy package.

**Key Data Structures:**

**1. Series (1D):**
```
+-------+-------+
| Index | Value |
+-------+-------+
|   0   | val_0 |
+-------+-------+
|   1   | val_1 |
+-------+-------+
|   2   | val_2 |
+-------+-------+
```

**2. DataFrame (2D):**
```
+-------+-------+-------+
| Index | Col A | Col B |
+-------+-------+-------+
|   0   |   1   |   2   |
+-------+-------+-------+
|   1   |   3   |   4   |
+-------+-------+-------+
|   2   |   5   |   6   |
+-------+-------+-------+
```

---

### **2.11.2 Pandas Series**

**Creating Series:**

```python
import pandas as pd
import numpy as np

# (i) From a list
data = pd.Series([0.25, 0.5, 0.75])

# With explicit index
data = pd.Series([0.25, 0.5, 0.75], index=['a', 'b', 'c'])
# data['b'] → 0.5

# (ii) From NumPy array
a_series = pd.Series(np.array([10, 20, 30]))

# (iii) From dictionary
d = {'one': 1, 'two': 2, 'three': 3}
s = pd.Series(d)
```

**Key Series Methods:**

```python
s.head(3)        # First 3 rows
s.tail(3)        # Last 3 rows
s.sort_values()  # Sort by values
s.sort_index()   # Sort by index
s.size           # Number of elements
s.is_unique      # True if all values unique
s.idxmax()       # Index of max value
```

---

### **2.11.3 Pandas DataFrame**

**Creating DataFrame:**

```python
import pandas as pd

raw_data = {
    'firstname': ['Rupa', 'Rakshita', 'Rupal'],
    'lastname': ['deri', 'Mani', 'deri'],
    'RNO': [12, 22, 12],
    'Testscore': [12, 22, 12]
}

df = pd.DataFrame(raw_data)

# Drop duplicates
df.drop_duplicates()  # Removes duplicate rows
```

**Example:**

```python
scottish_hills = {
    'Ben Nevis': (1345, 5),
    'Ben Macdui': (1405, 5)
}

dataframe = pd.DataFrame(scottish_hills)
```

---
---

### **2.12 PANDAS: SERIES AND DATAFRAME**

#### **2.12.1 Pandas Series**

**Definition:** A one-dimensional labeled array capable of holding any data type.

**Creating Series:**

```python
import pandas as pd
import numpy as np

# Method 1: From a list
data = pd.Series([0.25, 0.5, 0.75, 1.0])

# Method 2: With explicit index
data = pd.Series([0.25, 0.5, 0.75, 1.0], 
                 index=['a', 'b', 'c', 'd'])
print(data['b'])  # → 0.5

# Method 3: From NumPy array
arr = np.array([10, 20, 30, 40])
series = pd.Series(arr)

# Method 4: From dictionary
dict_data = {'one': 1, 'two': 2, 'three': 3}
s = pd.Series(dict_data)
```

**Key Series Methods:**

```python
s.head(3)           # First 3 rows
s.tail(3)           # Last 3 rows
s.sort_values()     # Sort by values
s.sort_index()      # Sort by index
s.size              # Number of elements
s.is_unique         # True if all values unique
s.idxmax()          # Index of maximum value
s.describe()        # Statistical summary
```

---

#### **2.12.2 Pandas DataFrame**

**Definition:** A two-dimensional labeled data structure with columns of potentially different types.

**Creating DataFrame:**

```python
import pandas as pd

# Method 1: From dictionary of lists
raw_data = {
    'firstname': ['Rupa', 'Rakshita', 'Rupal'],
    'lastname': ['deri', 'Mani', 'deri'],
    'RNO': [12, 22, 12],
    'Testscore': [12, 22, 12]
}

df = pd.DataFrame(raw_data)

# Method 2: Drop duplicates
df.drop_duplicates()  # Removes duplicate rows

# Method 3: From nested dictionary
scottish_hills = {
    'Ben Nevis': (1345, 5),
    'Ben Macdui': (1405, 5)
}
dataframe = pd.DataFrame(scottish_hills)
```

**ASCII Visualization:**

```
DataFrame Structure:
+----+-----------+----------+-----+-----------+
|    | firstname | lastname | RNO | Testscore |
+----+-----------+----------+-----+-----------+
| 0  | Rupa      | deri     | 12  | 12        |
| 1  | Rakshita  | Mani     | 22  | 22        |
| 2  | Rupal     | deri     | 12  | 12        |
+----+-----------+----------+-----+-----------+
```

---

## 2.13 MISSING DATA HANDLING

### **2.13.1 Understanding Missing Data**

**Definition:** In Python/Pandas, missing values are marked as **NaN (Not a Number)**. They are ignored from standard mathematical operations.

**Why Handle Missing Data?**
- Missing data can distort statistical analysis
- Most machine learning algorithms cannot handle NaN values
- Can lead to biased results

---

### **2.13.2 Detecting Missing Data**

```python
import pandas as pd
import numpy as np

df = pd.DataFrame({
    'A': [1, 2, np.nan, 4],
    'B': [5, np.nan, 7, 8],
    'C': [10, 11, 12, 13]
})

# Check for missing values (returns boolean DataFrame)
df.isnull()
# Output:
#        A      B      C
# 0  False  False  False
# 1  False   True  False
# 2   True  False  False
# 3  False  False  False

# Count missing values per column
df.isnull().sum()
# Output:
# A    1
# B    1
# C    0

# Total missing values
df.isnull().sum().sum()  # → 2
```

---

### **2.13.3 Methods to Handle Missing Data**

#### **Method 1: Fill with Constant Value**

```python
df.fillna(0)  # Replace all NaN with 0
```

#### **Method 2: Forward Fill (ffill)**

```python
# Fill NaN with the previous valid value
df.fillna(method='ffill')
# Or equivalently:
df.fillna(method='pad')
```

**Example:**
```
Before:          After ffill:
[1, NaN, 3]  →   [1, 1, 3]
[NaN, 5, 6]  →   [NaN, 5, 6]  (no previous value in row)
```

#### **Method 3: Backward Fill (bfill)**

```python
# Fill NaN with the next valid value
df.fillna(method='bfill')
# Or equivalently:
df.fillna(method='backfill')
```

**Example:**
```
Before:          After bfill:
[1, NaN, 3]  →   [1, 3, 3]
[NaN, 5, 6]  →   [5, 5, 6]
```

#### **Method 4: Fill with Column Mean/Median**

```python
# Fill with column mean (statistically neutral)
df.fillna(df.mean())

# Fill with column median (better for skewed data)
df.fillna(df.median())
```

**Example:**
```python
df = pd.DataFrame({
    'Age': [25, 30, np.nan, 40, 35]
})

mean_age = df['Age'].mean()  # → 32.5
df_filled = df.fillna(mean_age)
# Result: [25, 30, 32.5, 40, 35]
```

#### **Method 5: Drop Rows with NaN**

```python
# Drop any row containing at least one NaN
df.dropna()

# Drop only if ALL values are NaN
df.dropna(how='all')

# Drop only if NaN in specific column
df.dropna(subset=['A'])
```

#### **Method 6: Drop Columns with NaN**

```python
# Drop columns containing any NaN
df.dropna(axis=1)
```

---

### **2.13.4 Choosing the Right Method**

| Method | When to Use | Example |
|--------|-------------|---------|
| **Forward Fill** | Time-series data | Stock prices, temperature readings |
| **Backward Fill** | Time-series with future data | Forecasting models |
| **Mean/Median Fill** | Numerical continuous data | Age, income, test scores |
| **Mode Fill** | Categorical data | Gender, city, product category |
| **Drop Rows** | When missing data is minimal (<5%) | Small datasets with few missing values |
| **Drop Columns** | When column has too many missing values (>50%) | Irrelevant features |

---

## 2.14 HIERARCHICAL INDEXING (MULTIINDEX)

### **2.14.1 Definition**

**Hierarchical Indexing (MultiIndex)** enables storing and manipulating higher-dimensional data (3D, 4D) within the standard 1D Series and 2D DataFrame structures.

**Purpose:**
- Organize data with multiple levels of categories
- Efficient grouping and aggregation
- Reshape data easily

---

### **2.14.2 MultiIndexed Series**

```python
import pandas as pd
import numpy as np

# Create MultiIndex Series
index = [
    ["Maths", "Maths", "Maths", "Science", "Science", "Science"],
    ["Test1", "Test2", "Test3", "Test1", "Test2", "Test3"]
]

data = pd.Series([85, 90, 88, 78, 82, 80], index=index)

print(data)
# Output:
# Maths    Test1    85
#          Test2    90
#          Test3    88
# Science  Test1    78
#          Test2    82
#          Test3    80
# dtype: int64

# Access data by outer index
print(data["Maths"])
# Output:
# Test1    85
# Test2    90
# Test3    88

# Access specific value
print(data["Maths"]["Test1"])  # → 85

# Partial slicing
print(data["Maths":"Science"])
```

**ASCII Visualization:**

```
MultiIndex Series Structure:
Outer Index (Subject) → Maths
                          ↓
Inner Index (Test)   → Test1  Test2  Test3
                          ↓      ↓      ↓
Values               →   85     90     88
```

---

### **2.14.3 Unstack and Stack Operations**

**Unstack:** Converts MultiIndex Series to 2D DataFrame

```python
# Unstack inner level to columns
df = data.unstack()

print(df)
# Output:
#         Test1  Test2  Test3
# Maths      85     90     88
# Science    78     82     80

# Stack back to Series
df.stack()
```

---

### **2.14.4 MultiIndexed DataFrame**

```python
# Create MultiIndex DataFrame
arrays = [
    ['A', 'A', 'B', 'B'],
    [1, 2, 1, 2]
]

index = pd.MultiIndex.from_arrays(arrays, names=['Class', 'Exam'])

df = pd.DataFrame({
    'Python': [85, 90, 78, 82],
    'DS': [88, 92, 80, 85],
    'CA': [90, 88, 85, 87]
}, index=index)

print(df)
# Output:
#            Python   DS   CA
# Class Exam                
# A     1        85   88   90
#       2        90   92   88
# B     1        78   80   85
#       2        82   85   87

# Sort by specific level
df.sort_index(level='Exam')

# Access by outer index
df.loc['A']

# Access by both indices
df.loc[('A', 1)]
```

---

### **2.14.5 Practical Use Case: Student Exam Scores**

**Scenario:** Store student exam scores across multiple subjects and multiple test attempts.

```python
# MultiIndex: Outer = Subject, Inner = Test Number
index = [
    ["Maths", "Maths", "Science", "Science", "English", "English"],
    ["Test1", "Test2", "Test1", "Test2", "Test1", "Test2"]
]

scores = pd.Series([85, 90, 78, 88, 92, 87], index=index)

# Access all Maths scores
print(scores["Maths"])

# Access specific test
print(scores["Science"]["Test2"])  # → 88

# Convert to DataFrame for easier viewing
scores_df = scores.unstack()
print(scores_df)
# Output:
#         Test1  Test2
# Maths      85     90
# Science    78     88
# English    92     87

# Calculate average per subject
subject_avg = scores.groupby(level=0).mean()
print(subject_avg)
# Output:
# English    89.5
# Maths      87.5
# Science    83.0
```

---

## 2.15 COMBINING DATASETS

### **2.15.1 Concatenation**

**Definition:** Joining two or more DataFrames along a particular axis (rows or columns).

```python
import pandas as pd

# Concatenate Series
ser1 = pd.Series(['A', 'B', 'C'], index=[1, 2, 3])
ser2 = pd.Series(['X', 'Y', 'Z'], index=[4, 5, 6])

result = pd.concat([ser1, ser2])
print(result)
# Output:
# 1    A
# 2    B
# 3    C
# 4    X
# 5    Y
# 6    Z

# Concatenate DataFrames (vertical - axis=0)
df1 = pd.DataFrame({'Language': ['Python', 'NumPy', 'Pandas']})
df2 = pd.DataFrame({'Language': ['C', 'C++', 'Java']})

df = pd.concat([df1, df2])
print(df)
# Output:
#   Language
# 0   Python
# 1    NumPy
# 2   Pandas
# 0        C
# 1      C++
# 2     Java

# Concatenate horizontally (axis=1)
df3 = pd.DataFrame({'A': [1, 2, 3]})
df4 = pd.DataFrame({'B': [4, 5, 6]})

df_horizontal = pd.concat([df3, df4], axis=1)
print(df_horizontal)
# Output:
#    A  B
# 0  1  4
# 1  2  5
# 2  3  6
```

---

### **2.15.2 Append Method**

```python
# Append one DataFrame to another
df1 = pd.DataFrame({'A': [1, 2], 'B': [3, 4]})
df2 = pd.DataFrame({'A': [5, 6], 'B': [7, 8]})

result = df1.append(df2, ignore_index=True)
print(result)
# Output:
#    A  B
# 0  1  3
# 1  2  4
# 2  5  7
# 3  6  8
```

---

### **2.15.3 Merge/Join (SQL-style Joins)**

**Types of Joins:**

```python
df1 = pd.DataFrame({
    'key': ['A', 'B', 'C'],
    'value1': [1, 2, 3]
})

df2 = pd.DataFrame({
    'key': ['A', 'B', 'D'],
    'value2': [4, 5, 6]
})

# Inner Join (only matching keys)
pd.merge(df1, df2, on='key', how='inner')
# Output:
#   key  value1  value2
# 0   A       1       4
# 1   B       2       5

# Left Join (all from left, matching from right)
pd.merge(df1, df2, on='key', how='left')
# Output:
#   key  value1  value2
# 0   A       1     4.0
# 1   B       2     5.0
# 2   C       3     NaN

# Right Join (all from right, matching from left)
pd.merge(df1, df2, on='key', how='right')
# Output:
#   key  value1  value2
# 0   A     1.0       4
# 1   B     2.0       5
# 2   D     NaN       6

# Outer Join (all records from both)
pd.merge(df1, df2, on='key', how='outer')
# Output:
#   key  value1  value2
# 0   A     1.0     4.0
# 1   B     2.0     5.0
# 2   C     3.0     NaN
# 3   D     NaN     6.0
```

---

## 2.16 AGGREGATION AND GROUPBY

### **2.16.1 GroupBy: Split-Apply-Combine Strategy**

**Definition:** GroupBy follows three steps:
1. **Split** the DataFrame into groups based on a key column
2. **Apply** an aggregation function to each group
3. **Combine** the results back into a new DataFrame

---

### **2.16.2 Basic GroupBy Operations**

```python
import pandas as pd

# Create sample DataFrame
df = pd.DataFrame({
    'Dept': ['Sales', 'IT', 'Sales', 'IT', 'HR', 'HR'],
    'Salary': [50000, 75000, 48000, 82000, 52000, 55000],
    'Experience': [3, 7, 2, 9, 4, 6]
})

print(df)
# Output:
#     Dept  Salary  Experience
# 0  Sales   50000           3
# 1     IT   75000           7
# 2  Sales   48000           2
# 3     IT   82000           9
# 4     HR   52000           4
# 5     HR   55000           6

# Group by Department
grouped = df.groupby('Dept')

# Sum of salaries per department
print(grouped['Salary'].sum())
# Output:
# Dept
# HR     107000
# IT     157000
# Sales   98000

# Mean salary per department
print(grouped['Salary'].mean())
# Output:
# Dept
# HR     53500.0
# IT     78500.0
# Sales  49000.0

# Multiple aggregations
print(grouped['Salary'].agg(['mean', 'min', 'max', 'count']))
# Output:
#         mean    min    max  count
# Dept                              
# HR   53500.0  52000  55000      2
# IT   78500.0  75000  82000      2
# Sales 49000.0  48000  50000      2

# Group by multiple columns
df2 = df.copy()
df2['Year'] = [2022, 2022, 2023, 2023, 2022, 2023]

print(df2.groupby(['Dept', 'Year'])['Salary'].sum())
# Output:
# Dept   Year
# HR     2022    52000
#        2023    55000
# IT     2022    75000
#        2023    82000
# Sales  2022    50000
#        2023    48000
```

---

### **2.16.3 Common Aggregation Functions**

| Function | Description | Example |
|----------|-------------|---------|
| `sum()` | Sum of values | `grouped['Salary'].sum()` |
| `mean()` | Average | `grouped['Salary'].mean()` |
| `min()` | Minimum value | `grouped['Salary'].min()` |
| `max()` | Maximum value | `grouped['Salary'].max()` |
| `count()` | Count non-null values | `grouped['Salary'].count()` |
| `std()` | Standard deviation | `grouped['Salary'].std()` |
| `var()` | Variance | `grouped['Salary'].var()` |
| `first()` | First value in group | `grouped['Salary'].first()` |
| `last()` | Last value in group | `grouped['Salary'].last()` |
| `size()` | Total count (including NaN) | `grouped.size()` |
| `describe()` | Statistical summary | `grouped['Salary'].describe()` |

---

### **2.16.4 Practical Example: Sales Analysis**

```python
# Sales data
sales_df = pd.DataFrame({
    'Region': ['North', 'North', 'South', 'South', 'East', 'East'],
    'Product': ['A', 'B', 'A', 'B', 'A', 'B'],
    'Sales': [100, 200, 150, 250, 180, 220]
})

# Total sales per region
region_sales = sales_df.groupby('Region')['Sales'].sum()
print(region_sales)
# Output:
# Region
# East     400
# North    300
# South    400

# Average sales per product
product_avg = sales_df.groupby('Product')['Sales'].mean()
print(product_avg)
# Output:
# Product
# A    143.33
# B    223.33

# Multiple aggregations
summary = sales_df.groupby('Region')['Sales'].agg(['sum', 'mean', 'count'])
print(summary)
```

---

## 2.17 PIVOT TABLES

### **2.17.1 Definition**

**Pivot Tables** reorganize and summarize data from a column-wise format into a 2D tabular format, making it easy to compare groups along two dimensions simultaneously.

**Purpose:**
- Cross-tabulation of data
- Multi-dimensional analysis
- Summarize large datasets

---

### **2.17.2 Basic Pivot Table**

```python
import pandas as pd

# Sales data
sales_data = pd.DataFrame({
    'Region': ['North', 'North', 'South', 'South', 'East', 'East'],
    'Product': ['A', 'B', 'A', 'B', 'A', 'B'],
    'Year': [2022, 2022, 2022, 2022, 2023, 2023],
    'Sales': [100, 200, 150, 250, 180, 220]
})

# Basic pivot table
pivot = pd.pivot_table(
    sales_data,
    index='Region',           # Rows
    columns='Product',        # Columns
    values='Sales',           # Values to aggregate
    aggfunc='sum'             # Aggregation function
)

print(pivot)
# Output:
# Product    A    B
# Region           
# East     180  220
# North    100  200
# South    150  250
```

---

### **2.17.3 Advanced Pivot Tables**

**Multi-index Pivot:**

```python
# Pivot with multiple index levels
pivot2 = pd.pivot_table(
    sales_data,
    index=['Region', 'Year'],
    values='Sales',
    aggfunc=['sum', 'mean']
)

print(pivot2)
# Output:
#              sum    mean
# Region Year             
# East   2023  400   200.0
# North  2022  300   150.0
# South  2022  400   200.0
```

**With Margins (Totals):**

```python
pivot_with_totals = pd.pivot_table(
    sales_data,
    index='Region',
    columns='Product',
    values='Sales',
    aggfunc='sum',
    margins=True,           # Add row/column totals
    margins_name='Total'    # Name for totals
)

print(pivot_with_totals)
# Output:
# Product    A    B  Total
# Region                 
# East     180  220    400
# North    100  200    300
# South    150  250    400
# Total    430  670   1100
```

---

### **2.17.4 Practical Example: Multi-Year Sales Analysis**

```python
# Create comprehensive dataset
data = {
    'Year': [2021, 2021, 2021, 2021, 2022, 2022, 2022, 2022, 2023, 2023, 2023, 2023],
    'Region': ['North', 'South', 'East', 'West'] * 3,
    'Category': ['Electronics', 'Clothing', 'Electronics', 'Clothing',
                 'Electronics', 'Clothing', 'Electronics', 'Clothing',
                 'Electronics', 'Clothing', 'Electronics', 'Clothing'],
    'Sales': [500, 300, 450, 200, 600, 350, 480, 220, 700, 400, 520, 280]
}

df = pd.DataFrame(data)

# Pivot: Region × Year
pivot_region_year = pd.pivot_table(
    df,
    index='Region',
    columns='Year',
    values='Sales',
    aggfunc='sum',
    margins=True,
    margins_name='Total'
)

print(pivot_region_year)
# Output:
# Year     2021  2022  2023  Total
# Region                         
# East      450   480   520   1450
# North     500   600   700   1800
# South     300   350   400   1050
# West      200   220   280    700
# Total    1450  1650  1900   5000

# Pivot: Region × Category
pivot_category = pd.pivot_table(
    df,
    index='Region',
    columns='Category',
    values='Sales',
    aggfunc='sum'
)

print(pivot_category)
# Output:
# Category  Clothing  Electronics
# Region                        
# East           300         1150
# North          400         1400
# South          350          700
# West           200          500

# Insights:
# - North region is highest revenue generator (₹1800)
# - Electronics outsells Clothing in all regions
# - Sales growing year-on-year (1450→1650→1900)
```

---

### **2.17.5 Comparison: GroupBy vs Pivot Table**

| Feature | GroupBy | Pivot Table |
|---------|---------|-------------|
| **Output Shape** | 1D Series or DataFrame | 2D cross-tabular format |
| **Best For** | Single-dimension grouping | Multi-dimensional comparison |
| **Flexibility** | More flexible aggregations | Easier visualization |
| **Syntax** | `df.groupby('col').agg()` | `pd.pivot_table()` |
| **Example** | Sales by Region | Sales by Region × Product |

---
---

### **2.18 COMPARISON: GROUPBY VS PIVOT TABLE**

| Feature | GroupBy | Pivot Table |
|---------|---------|-------------|
| **Output Shape** | 1D Series or DataFrame | 2D cross-tabular format |
| **Best For** | Single-dimension grouping | Multi-dimensional comparison |
| **Flexibility** | More flexible aggregations | Easier visualization |
| **Syntax** | `df.groupby('col').agg()` | `pd.pivot_table()` |
| **Example** | Sales by Region | Sales by Region × Product |

---

## 📊 UNIT 3: DATA WRANGLING AND DATA VISUALIZATION

### **3.1 DATA INDEXING AND SELECTION IN PANDAS**

#### **3.1.1 Label-Based Indexing (.loc)**

Accesses data by **row and column labels**.

```python
import pandas as pd

df = pd.DataFrame({
    'Name': ['Alice', 'Bob', 'Carol', 'David'],
    'Dept': ['Sales', 'IT', 'Sales', 'IT'],
    'Salary': [50000, 75000, 48000, 82000],
    'Experience': [3, 7, 2, 9]
})

# Select single row by label
print(df.loc[0])
# Output:
# Name          Alice
# Dept          Sales
# Salary        50000
# Experience        3

# Select range of rows and columns
print(df.loc[0:2, 'Name':'Salary'])
# Output:
#     Name   Dept  Salary
# 0  Alice  Sales   50000
# 1    Bob     IT   75000
# 2  Carol  Sales   48000

# Select specific rows and columns
print(df.loc[[0, 2], ['Name', 'Salary']])
# Output:
#     Name  Salary
# 0  Alice   50000
# 2  Carol   48000
```

---

#### **3.1.2 Position-Based Indexing (.iloc)**

Accesses data by **integer position** (0-based indexing).

```python
# Select first row by position
print(df.iloc[0])
# Output: Same as df.loc[0]

# Select first 2 rows, first 3 columns
print(df.iloc[0:2, 0:3])
# Output:
#     Name   Dept  Salary
# 0  Alice  Sales   50000
# 1    Bob     IT   75000

# Select last row
print(df.iloc[-1])
# Output: David's record

# Select specific positions
print(df.iloc[[0, 3], [1, 3]])
# Output:
#     Dept  Experience
# 0  Sales           3
# 3     IT           9
```

---

#### **3.1.3 Boolean Indexing**

Filter rows based on **conditions**.

```python
# Filter: Sales department with salary > 45000
result = df[(df['Dept'] == 'Sales') & (df['Salary'] > 45000)]
print(result)
# Output:
#     Name   Dept  Salary  Experience
# 0  Alice  Sales   50000           3

# Filter: Experience > 5 OR Salary > 70000
result2 = df[(df['Experience'] > 5) | (df['Salary'] > 70000)]
print(result2)
# Output:
#   Name Dept  Salary  Experience
# 1  Bob   IT   75000           7
# 3 David   IT   82000           9

# Filter: Salary between 50000 and 80000
result3 = df[df['Salary'].between(50000, 80000)]
print(result3)
# Output:
#     Name   Dept  Salary  Experience
# 0  Alice  Sales   50000           3
# 1    Bob     IT   75000           7
```

---

### **3.2 OPERATING ON DATA**

#### **3.2.1 Vectorized Operations**

Perform operations on entire columns without loops.

```python
# Add 10% raise to all salaries
df['New_Salary'] = df['Salary'] * 1.10
print(df[['Salary', 'New_Salary']])
# Output:
#    Salary  New_Salary
# 0   50000     55000.0
# 1   75000     82500.0
# 2   48000     52800.0
# 3   82000     90200.0

# Calculate salary per year of experience
df['Salary_Per_Year'] = df['Salary'] / df['Experience']
print(df[['Experience', 'Salary', 'Salary_Per_Year']])
# Output:
#    Experience  Salary  Salary_Per_Year
# 0           3   50000        16666.67
# 1           7   75000        10714.29
# 2           2   48000        24000.00
# 3           9   82000         9111.11

# String operations
df['Name_Upper'] = df['Name'].str.upper()
print(df['Name_Upper'])
# Output:
# 0    ALICE
# 1      BOB
# 2    CAROL
# 3    DAVID
```

---

#### **3.2.2 Applying Functions**

```python
# Apply custom function to a column
def salary_category(salary):
    if salary < 50000:
        return 'Low'
    elif salary < 75000:
        return 'Medium'
    else:
        return 'High'

df['Category'] = df['Salary'].apply(salary_category)
print(df[['Salary', 'Category']])
# Output:
#    Salary  Category
# 0   50000    Medium
# 1   75000      High
# 2   48000       Low
# 3   82000      High

# Apply to multiple columns
df[['Salary', 'Experience']].apply(np.mean)
# Output:
# Salary        63750.0
# Experience        5.25
```

---

### **3.3 MISSING DATA HANDLING**

#### **3.3.1 Understanding Missing Data**

Missing values in Pandas are represented as **NaN (Not a Number)** or **None**.

```python
import numpy as np

df = pd.DataFrame({
    'A': [1, 2, np.nan, 4],
    'B': [5, np.nan, 7, 8],
    'C': [10, 11, 12, 13]
})

print(df)
# Output:
#      A    B   C
# 0  1.0  5.0  10
# 1  2.0  NaN  11
# 2  NaN  7.0  12
# 3  4.0  8.0  13
```

---

#### **3.3.2 Detecting Missing Values**

```python
# Check for missing values (returns boolean DataFrame)
print(df.isnull())
# Output:
#        A      B      C
# 0  False  False  False
# 1  False   True  False
# 2   True  False  False
# 3  False  False  False

# Count missing values per column
print(df.isnull().sum())
# Output:
# A    1
# B    1
# C    0

# Total missing values
print(df.isnull().sum().sum())  # → 2

# Check for non-null values
print(df.notnull())
```

---

#### **3.3.3 Methods to Handle Missing Data**

**Method 1: Fill with Constant Value**

```python
df_filled = df.fillna(0)
print(df_filled)
# Output:
#      A    B   C
# 0  1.0  5.0  10
# 1  2.0  0.0  11
# 2  0.0  7.0  12
# 3  4.0  8.0  13
```

**Method 2: Forward Fill (ffill)**

```python
# Fill NaN with the previous valid value
df_ffill = df.fillna(method='ffill')
print(df_ffill)
# Output:
#      A    B   C
# 0  1.0  5.0  10
# 1  2.0  5.0  11
# 2  2.0  7.0  12
# 3  4.0  8.0  13
```

**Method 3: Backward Fill (bfill)**

```python
# Fill NaN with the next valid value
df_bfill = df.fillna(method='bfill')
print(df_bfill)
# Output:
#      A    B   C
# 0  1.0  5.0  10
# 1  2.0  7.0  11
# 2  4.0  7.0  12
# 3  4.0  8.0  13
```

**Method 4: Fill with Column Mean/Median**

```python
# Fill with column mean
df_mean = df.fillna(df.mean())
print(df_mean)
# Output:
#      A    B   C
# 0  1.0  5.0  10
# 1  2.0  6.67 11
# 2  2.33 7.0  12
# 3  4.0  8.0  13

# Fill with column median
df_median = df.fillna(df.median())
print(df_median)
```

**Method 5: Drop Rows with NaN**

```python
# Drop any row containing at least one NaN
df_dropped = df.dropna()
print(df_dropped)
# Output:
#      A    B   C
# 0  1.0  5.0  10
# 3  4.0  8.0  13

# Drop only if ALL values are NaN
df.dropna(how='all')

# Drop only if NaN in specific column
df.dropna(subset=['A'])
```

**Method 6: Drop Columns with NaN**

```python
# Drop columns containing any NaN
df.dropna(axis=1)
```

---

#### **3.3.4 Choosing the Right Method**

| Method | When to Use | Example |
|--------|-------------|---------|
| **Forward Fill** | Time-series data | Stock prices, temperature readings |
| **Backward Fill** | Time-series with future data | Forecasting models |
| **Mean/Median Fill** | Numerical continuous data | Age, income, test scores |
| **Mode Fill** | Categorical data | Gender, city, product category |
| **Drop Rows** | When missing data is minimal (<5%) | Small datasets with few missing values |
| **Drop Columns** | When column has too many missing values (>50%) | Irrelevant features |

---

### **3.4 HIERARCHICAL INDEXING (MULTIINDEX)**

#### **3.4.1 MultiIndexed Series**

```python
import pandas as pd
import numpy as np

# Create MultiIndex Series
index = [
    ["Maths", "Maths", "Maths", "Science", "Science", "Science"],
    ["Test1", "Test2", "Test3", "Test1", "Test2", "Test3"]
]

data = pd.Series([85, 90, 88, 78, 82, 80], index=index)

print(data)
# Output:
# Maths    Test1    85
#          Test2    90
#          Test3    88
# Science  Test1    78
#          Test2    82
#          Test3    80

# Access data by outer index
print(data["Maths"])
# Output:
# Test1    85
# Test2    90
# Test3    88

# Access specific value
print(data["Maths"]["Test1"])  # → 85

# Partial slicing
print(data["Maths":"Science"])
```

**ASCII Visualization:**

```
MultiIndex Series Structure:
Outer Index (Subject) → Maths
                          ↓
Inner Index (Test)   → Test1  Test2  Test3
                          ↓      ↓      ↓
Values               →   85     90     88
```

---

#### **3.4.2 Unstack and Stack Operations**

**Unstack:** Converts MultiIndex Series to 2D DataFrame

```python
# Unstack inner level to columns
df = data.unstack()

print(df)
# Output:
#         Test1  Test2  Test3
# Maths      85     90     88
# Science    78     82     80

# Stack back to Series
df.stack()
```

---

#### **3.4.3 MultiIndexed DataFrame**

```python
# Create MultiIndex DataFrame
arrays = [
    ['A', 'A', 'B', 'B'],
    [1, 2, 1, 2]
]

index = pd.MultiIndex.from_arrays(arrays, names=['Class', 'Exam'])

df = pd.DataFrame({
    'Python': [85, 90, 78, 82],
    'DS': [88, 92, 80, 85],
    'CA': [90, 88, 85, 87]
}, index=index)

print(df)
# Output:
#            Python   DS   CA
# Class Exam                
# A     1        85   88   90
#       2        90   92   88
# B     1        78   80   85
#       2        82   85   87

# Sort by specific level
df.sort_index(level='Exam')

# Access by outer index
df.loc['A']

# Access by both indices
df.loc[('A', 1)]
```

---

### **3.5 COMBINING DATASETS**

#### **3.5.1 Concatenation**

```python
import pandas as pd

# Concatenate Series
ser1 = pd.Series(['A', 'B', 'C'], index=[1, 2, 3])
ser2 = pd.Series(['X', 'Y', 'Z'], index=[4, 5, 6])

result = pd.concat([ser1, ser2])
print(result)
# Output:
# 1    A
# 2    B
# 3    C
# 4    X
# 5    Y
# 6    Z

# Concatenate DataFrames (vertical - axis=0)
df1 = pd.DataFrame({'Language': ['Python', 'NumPy', 'Pandas']})
df2 = pd.DataFrame({'Language': ['C', 'C++', 'Java']})

df = pd.concat([df1, df2])
print(df)
# Output:
#   Language
# 0   Python
# 1    NumPy
# 2   Pandas
# 0        C
# 1      C++
# 2     Java

# Concatenate horizontally (axis=1)
df3 = pd.DataFrame({'A': [1, 2, 3]})
df4 = pd.DataFrame({'B': [4, 5, 6]})

df_horizontal = pd.concat([df3, df4], axis=1)
print(df_horizontal)
# Output:
#    A  B
# 0  1  4
# 1  2  5
# 2  3  6
```

---

#### **3.5.2 Append Method**

```python
# Append one DataFrame to another
df1 = pd.DataFrame({'A': [1, 2], 'B': [3, 4]})
df2 = pd.DataFrame({'A': [5, 6], 'B': [7, 8]})

result = df1.append(df2, ignore_index=True)
print(result)
# Output:
#    A  B
# 0  1  3
# 1  2  4
# 2  5  7
# 3  6  8
```

---

#### **3.5.3 Merge/Join (SQL-style Joins)**

```python
df1 = pd.DataFrame({
    'key': ['A', 'B', 'C'],
    'value1': [1, 2, 3]
})

df2 = pd.DataFrame({
    'key': ['A', 'B', 'D'],
    'value2': [4, 5, 6]
})

# Inner Join (only matching keys)
pd.merge(df1, df2, on='key', how='inner')
# Output:
#   key  value1  value2
# 0   A       1       4
# 1   B       2       5

# Left Join (all from left, matching from right)
pd.merge(df1, df2, on='key', how='left')
# Output:
#   key  value1  value2
# 0   A       1     4.0
# 1   B       2     5.0
# 2   C       3     NaN

# Right Join (all from right, matching from left)
pd.merge(df1, df2, on='key', how='right')
# Output:
#   key  value1  value2
# 0   A     1.0       4
# 1   B     2.0       5
# 2   D     NaN       6

# Outer Join (all records from both)
pd.merge(df1, df2, on='key', how='outer')
# Output:
#   key  value1  value2
# 0   A     1.0     4.0
# 1   B     2.0     5.0
# 2   C     3.0     NaN
# 3   D     NaN     6.0
```

---

### **3.6 AGGREGATION AND GROUPBY**

#### **3.6.1 GroupBy: Split-Apply-Combine Strategy**

**Definition:** GroupBy follows three steps:
1. **Split** the DataFrame into groups based on a key column
2. **Apply** an aggregation function to each group
3. **Combine** the results back into a new DataFrame

---

#### **3.6.2 Basic GroupBy Operations**

```python
import pandas as pd

# Create sample DataFrame
df = pd.DataFrame({
    'Dept': ['Sales', 'IT', 'Sales', 'IT', 'HR', 'HR'],
    'Salary': [50000, 75000, 48000, 82000, 52000, 55000],
    'Experience': [3, 7, 2, 9, 4, 6]
})

print(df)
# Output:
#     Dept  Salary  Experience
# 0  Sales   50000           3
# 1     IT   75000           7
# 2  Sales   48000           2
# 3     IT   82000           9
# 4     HR   52000           4
# 5     HR   55000           6

# Group by Department
grouped = df.groupby('Dept')

# Sum of salaries per department
print(grouped['Salary'].sum())
# Output:
# Dept
# HR     107000
# IT     157000
# Sales   98000

# Mean salary per department
print(grouped['Salary'].mean())
# Output:
# Dept
# HR     53500.0
# IT     78500.0
# Sales  49000.0

# Multiple aggregations
print(grouped['Salary'].agg(['mean', 'min', 'max', 'count']))
# Output:
#         mean    min    max  count
# Dept                              
# HR   53500.0  52000  55000      2
# IT   78500.0  75000  82000      2
# Sales 49000.0  48000  50000      2

# Group by multiple columns
df2 = df.copy()
df2['Year'] = [2022, 2022, 2023, 2023, 2022, 2023]

print(df2.groupby(['Dept', 'Year'])['Salary'].sum())
# Output:
# Dept   Year
# HR     2022    52000
#        2023    55000
# IT     2022    75000
#        2023    82000
# Sales  2022    50000
#        2023    48000
```

---

#### **3.6.3 Common Aggregation Functions**

| Function | Description | Example |
|----------|-------------|---------|
| `sum()` | Sum of values | `grouped['Salary'].sum()` |
| `mean()` | Average | `grouped['Salary'].mean()` |
| `min()` | Minimum value | `grouped['Salary'].min()` |
| `max()` | Maximum value | `grouped['Salary'].max()` |
| `count()` | Count non-null values | `grouped['Salary'].count()` |
| `std()` | Standard deviation | `grouped['Salary'].std()` |
| `var()` | Variance | `grouped['Salary'].var()` |
| `first()` | First value in group | `grouped['Salary'].first()` |
| `last()` | Last value in group | `grouped['Salary'].last()` |
| `size()` | Total count (including NaN) | `grouped.size()` |
| `describe()` | Statistical summary | `grouped['Salary'].describe()` |

---

### **3.7 PIVOT TABLES**

#### **3.7.1 Definition**

**Pivot Tables** reorganize and summarize data from a column-wise format into a 2D tabular format, making it easy to compare groups along two dimensions simultaneously.

**Purpose:**
- Cross-tabulation of data
- Multi-dimensional analysis
- Summarize large datasets

---

#### **3.7.2 Basic Pivot Table**

```python
import pandas as pd

# Sales data
sales_data = pd.DataFrame({
    'Region': ['North', 'North', 'South', 'South', 'East', 'East'],
    'Product': ['A', 'B', 'A', 'B', 'A', 'B'],
    'Year': [2022, 2022, 2022, 2022, 2023, 2023],
    'Sales': [100, 200, 150, 250, 180, 220]
})

# Basic pivot table
pivot = pd.pivot_table(
    sales_data,
    index='Region',           # Rows
    columns='Product',        # Columns
    values='Sales',           # Values to aggregate
    aggfunc='sum'             # Aggregation function
)

print(pivot)
# Output:
# Product    A    B
# Region           
# East     180  220
# North    100  200
# South    150  250
```

---

#### **3.7.3 Advanced Pivot Tables**

**Multi-index Pivot:**

```python
# Pivot with multiple index levels
pivot2 = pd.pivot_table(
    sales_data,
    index=['Region', 'Year'],
    values='Sales',
    aggfunc=['sum', 'mean']
)

print(pivot2)
# Output:
#              sum    mean
# Region Year             
# East   2023  400   200.0
# North  2022  300   150.0
# South  2022  400   200.0
```

**With Margins (Totals):**

```python
pivot_with_totals = pd.pivot_table(
    sales_data,
    index='Region',
    columns='Product',
    values='Sales',
    aggfunc='sum',
    margins=True,           # Add row/column totals
    margins_name='Total'    # Name for totals
)

print(pivot_with_totals)
# Output:
# Product    A    B  Total
# Region                 
# East     180  220    400
# North    100  200    300
# South    150  250    400
# Total    430  670   1100
```

---

#### **3.7.4 Practical Example: Multi-Year Sales Analysis**

```python
# Create comprehensive dataset
data = {
    'Year': [2021, 2021, 2021, 2021, 2022, 2022, 2022, 2022, 2023, 2023, 2023, 2023],
    'Region': ['North', 'South', 'East', 'West'] * 3,
    'Category': ['Electronics', 'Clothing', 'Electronics', 'Clothing',
                 'Electronics', 'Clothing', 'Electronics', 'Clothing',
                 'Electronics', 'Clothing', 'Electronics', 'Clothing'],
    'Sales': [500, 300, 450, 200, 600, 350, 480, 220, 700, 400, 520, 280]
}

df = pd.DataFrame(data)

# Pivot: Region × Year
pivot_region_year = pd.pivot_table(
    df,
    index='Region',
    columns='Year',
    values='Sales',
    aggfunc='sum',
    margins=True,
    margins_name='Total'
)

print(pivot_region_year)
# Output:
# Year     2021  2022  2023  Total
# Region                         
# East      450   480   520   1450
# North     500   600   700   1800
# South     300   350   400   1050
# West      200   220   280    700
# Total    1450  1650  1900   5000

# Pivot: Region × Category
pivot_category = pd.pivot_table(
    df,
    index='Region',
    columns='Category',
    values='Sales',
    aggfunc='sum'
)

print(pivot_category)
# Output:
# Category  Clothing  Electronics
# Region                        
# East           300         1150
# North          400         1400
# South          350          700
# West           200          500

# Insights:
# - North region is highest revenue generator (₹1800)
# - Electronics outsells Clothing in all regions
# - Sales growing year-on-year (1450→1650→1900)
```

---

### **3.8 IMPORTING MATPLOTLIB**

#### **3.8.1 Basic Line Plots**

```python
import matplotlib.pyplot as plt
import numpy as np

# Basic plot
plt.plot([1, 2, 3], [5, 7, 4])
plt.show()

# Line Plots - multiple lines
x = np.linspace(-1, 1, 50)
y1 = 2 * x + 1          # linear line
y2 = 2 * x**2 + 1       # curved (non-linear)

plt.figure(num=3, figsize=(8, 5))
plt.plot(x, y2)
plt.plot(x, y1, linewidth=1.0, linestyle='--')
plt.show()

# With axis labels
x = np.array([1, 2, 3, 4])
y = x**2

plt.plot(x, y)
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title('Line Graph')
plt.show()

# Multi-point line: (1,3)→(2,4)→(6,1)→(8,10)
x = np.array([1, 2, 6, 8])
y = np.array([3, 4, 1, 10])
plt.plot(x, y)
plt.show()
```

---

#### **3.8.2 Saving Figures**

```python
# Saving Figures
plt.savefig('myplot.png', 
            dpi=300,              # High resolution
            transparent=True,     # Transparent background
            bbox_inches='tight')  # Remove extra whitespace
```

---

#### **3.8.3 Axes, Ticks, and Grids**

```python
# Axes, Ticks, Grids
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # [left, bottom, width, height]
axes.plot(x, y, 'r')  # 'r' = red color
axes.set_xlabel('x')
axes.set_ylabel('y')
axes.set_title('title')

plt.grid(True)  # enable grid
```

---

### **3.9 SCATTER PLOTS**

#### **3.9.1 Basic Scatter Plot**

```python
import matplotlib.pyplot as plt

# Basic scatter
x = [2, 3, 7, 29, 8]
y = [4, 7, 55, 43, 4]
plt.scatter(x, y)
plt.show()

# Advanced scatter - two datasets
x1 = [1, 2, 3, 4, 5]
y1 = [5, 16, 34, 56, 42]
x2 = [6, 7, 8, 9, 10]
y2 = [11, 12, 13, 14, 15]

plt.title("Prices over 10 years")
plt.scatter(x1, y1, color='darkblue', marker='x')
plt.scatter(x2, y2, color='red', marker='v')
plt.xlabel("Time")
plt.ylabel("Price")
plt.grid(True)
plt.legend()
plt.show()
```

---

#### **3.9.2 Common Markers**

| Marker | Description | Marker | Description |
|--------|-------------|--------|-------------|
| `^` | Triangle up | `v` | Triangle down |
| `<` | Triangle left | `>` | Triangle right |
| `*` | Star | `+` | Plus |
| `o` | Circle | `s` | Square |
| `S` | Filled square | `h` | Hexagon |

---

#### **3.9.3 Labels, Annotations, and Legends**

```python
# Annotations
plt.annotate("Maximum", 
             xy=(3, 3),           # Point to annotate
             xytext=(3, 1.8),     # Text position
             arrowprops=dict(facecolor='green'))

# Legend on bottom
ax.legend(loc='upper center', 
          bbox_to_anchor=(0.5, -0.5), 
          shadow=True, 
          ncol=2)
```

---

### **3.10 VISUALIZING ERRORS**

#### **3.10.1 Error Bars**

```python
import matplotlib.pyplot as plt
import numpy as np

# Basic error bar
x = np.arange(1, 8)
y = np.array([20, 10, 45, 32, 38, 21, 29])
yerror = y * 0.10  # 10% error

# fmt: o=filled circle, b=blue. ecolor= error bar color(k=black)
plt.errorbar(x, y, yerr=yerror, 
             linestyle='None', 
             fmt='ob', 
             ecolor='k', 
             capsize=3)
plt.show()

# Advanced limits
plt.errorbar(x, y + 8, yerr=yerror)
plt.errorbar(x, y + 6, yerr=yerror, uplims=True)
plt.errorbar(x, y + 4, yerr=yerror, uplims=True, lolims=True)
```

---

### **3.11 DENSITY AND CONTOUR PLOTS**

#### **3.11.1 Contour Plots**

```python
import numpy as np
import matplotlib.pyplot as plt

def f(x, y):
    return np.sin(x)**10 + np.cos(10 + y * x) * np.cos(x)

x = np.linspace(0, 5, 50)
y = np.linspace(0, 5, 40)
X, Y = np.meshgrid(x, y)  # creates 2D coordinate arrays
Z = f(X, Y)

plt.style.use('seaborn-white')
plt.contour(X, Y, Z, colors='black')  # line contour
plt.contourf(X, Y, Z)  # filled contour
plt.show()
```

**ASCII Visualization:**

```
Contour Plot (Topographic Map):
        Y
        ↑
    ┌───┼───┐
    │ ╱ ╱ ╱ │  ← Contour lines
    │╱ ╱ ╱ ╱│     (equal values)
    │╲ ╲ ╲ ╲│
    │ ╲ ╲ ╲ │
    └───┼───┘
        → X
```

---

#### **3.11.2 Histograms**

```python
x = 40 + np.random.randn(50000)

# bins=20, range=(-5, 50), align='mid', histtype='stepfilled'
plt.hist(x, 20, range=(-5, 50), 
         histtype='stepfilled', 
         align='mid', 
         color='r', 
         label='Test')
plt.legend()
plt.title('Histogram')
plt.show()
```

---

### **3.12 SUBPLOTS**

#### **3.12.1 Multiple Subplots**

```python
# Multiple subplots (grid arrangement)
for i in range(1, 7):
    plt.subplot(2, 3, i)  # 2 rows, 3 cols, position i
    plt.xticks([])
    plt.yticks([])
    plt.plot(x, y)
plt.show()
```

**ASCII Visualization:**

```
Subplot Grid (2 rows × 3 cols):
┌─────┬─────┬─────┐
│  1  │  2  │  3  │
├──────────┼─────┤
│  4  │  5  │  6  │
└─────┴─────┴─────┘
```

---

### **3.13 TEXT AND ANNOTATION**

```python
# Customization - Major/Minor Ticks
ax.minorticks_on()
ax.grid(which='both')  # show both major and minor grid

# Text annotation
plt.text(x, y, 'Sample Text', 
         fontsize=12, 
         color='red',
         bbox=dict(facecolor='white', alpha=0.5))
```

---

### **3.14 THREE-DIMENSIONAL PLOTTING**

```python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')  # 3D axes

z = np.linspace(0, 1, 10)
x = z * np.sin(20 * z)
y = z * np.cos(20 * z)
c = x + y

ax.scatter(x, y, z, c=c)  # colour by x+y
plt.show()
```

**ASCII Visualization:**

```
3D Scatter Plot:
           Z
           ↑
          /|\
         / | \
        /  |  \
       /   |   \
      /    ●    \
     /   ●   ●   \
    / ●         ● \
   └───────────────→ Y
  /
 /
X
```

---

### **3.15 GEOGRAPHIC DATA WITH BASEMAP**

**Basemap Methods:**

| Method | Purpose |
|--------|---------|
| `contour()` / `contourf()` | Draw contour lines / filled contours |
| `imshow()` / `pcolor()` | Draw image / pseudo-color plot |
| `plot()` / `scatter()` | Draw lines / points |
| `quiver()` / `barbs()` | Draw vectors / wind barbs |

---

### **3.16 VISUALIZATION WITH SEABORN**

#### **3.16.1 Matplotlib vs Seaborn**

| Feature | Matplotlib | Seaborn |
|---------|------------|---------|
| **Use Case** | Base library. Works with arrays. | Enhanced. Uses Matplotlib+NumPy+Pandas. |
| **Syntax** | Comparatively lengthy & complex | Comparatively simple & elegant |
| **Multiple Figs** | Can open multiple at a time | Automatically manages multiple figures |
| **Flexibility** | Highly customizable & powerful | Provides default production themes |

---

#### **3.16.2 Seaborn Examples**

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_csv("WorldHappiness2016.csv")

# Scatter plot
sns.scatterplot(data=df, x="Economy", y="Happiness")

# Line plot
sns.lineplot(data=df, x='x', y='y')

# Bar plot
sns.barplot(data=df, x='cat', y='val')

# Histogram
sns.histplot(data=df, x='col', bins=20)

# Box plot
sns.boxplot(data=df, x='cat', y='val')

# Heatmap (correlation matrix)
sns.heatmap(df.corr(), annot=True)

# Pair plot (all pairwise relationships)
sns.pairplot(df)

# Violin plot
sns.violinplot(data=df, x='cat', y='val')

# Set style
sns.set_style("whitegrid")  # darkgrid, dark, white, ticks
```

---

## 📝 QUICK REVISION CHEATSHEET

### **NumPy**

```python
np.array([1,2,3])              # Create 1D array
np.zeros((3,4))                # Array of zeros (3x4)
np.ones((2,3))                 # Array of ones
np.arange(0,10,2)              # 0, 2, 4, 6, 8
np.linspace(0,1,5)             # 5 evenly spaced in [0,1]
np.random.randn(3,3)           # 3x3 random normal array
arr.reshape(2,3)               # Reshape to 2 rows 3 cols
arr.T                          # Transpose
np.concatenate([a,b])          # Join arrays
np.split(arr,[3,5])            # Split at indices 3 & 5
np.sum/np.min/np.max/np.mean/np.std  # Aggregations
```

---

### **Pandas**

```python
pd.Series([1,2,3])             # Create Series
pd.DataFrame(data)             # Create DataFrame
df.head(5)/df.tail(5)          # First 5 rows/Last 5 rows
df.info()                      # Column types & NaN count
df.describe()                  # Statistical summary
df.isnull().sum()              # Count missing per column
df.fillna(0)                   # Fill NaN with 0
df.dropna()                    # Drop rows with NaN
df.drop_duplicates()           # Remove duplicate rows
df.groupby("col").sum()        # GroupBy and sum
pd.concat([df1,df2])           # Concatenate DataFrames
pd.merge(df1,df2,on='key')     # Merge (SQL join)
pd.pivot_table(df, index=["a"], aggfunc="mean")  # Pivot table
```

---

### **Matplotlib**

```python
import matplotlib.pyplot as plt

plt.plot(x,y)                  # Line plot
plt.scatter(x,y)               # Scatter plot
plt.hist(x,bins=20)            # Histogram
plt.errorbar(x, y, yerr=e)     # Error bar
plt.contour(X,Y,Z)             # Contour lines
plt.contourf(X,Y,Z)            # Filled contour
plt.xlabel("x"); plt.ylabel("y")  # Axis labels
plt.title("Title")             # Set title
plt.legend()                   # Show legend
plt.grid(True)                 # Enable grid
plt.show()                     # Display plot
plt.savefig("out.png",dpi=300) # Save figure
plt.subplot(rows,cols,pos)     # Create subplot
plt.annotate("text", xy, xytext, arrowprops)  # Annotation
from mpl_toolkits.mplot3d import Axes3D  # 3D support
```

---

### **Seaborn**

```python
import seaborn as sns

sns.scatterplot(data=df, x="x", y="y")  # Scatter
sns.lineplot(data=df, x="x", y="y")     # Line
sns.histplot(data=df, x="col")          # Histogram
sns.boxplot(data=df, x="cat", y="val")  # Box plot
sns.heatmap(df.corr(), annot=True)      # Heatmap
sns.pairplot(df)                        # Pair plot
sns.set_style("whitegrid")              # Set style
```

---



