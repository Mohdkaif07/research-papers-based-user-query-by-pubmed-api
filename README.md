# 🔍 Research Paper Fetching using PubMed API

This project is a Python-based CLI tool that allows users to search for research papers from the **PubMed database** based on a user-specified query. It filters the results to include only papers with at least one author affiliated with a **pharmaceutical** or **biotech** company and saves the output to a CSV file.

---

## 🔑 **Key Features**
♦️ **User-Specified Query** – Search for papers using any keyword or topic.  
♦️ **PubMed API Integration** – Fetch research papers directly from PubMed using `Entrez`.  
♦️ **Author Affiliation Filtering** – Filters papers based on pharma/biotech-related affiliations.  
♦️ **Structured CSV Output** – Saves results to a CSV file with key details.  
♦️ **Command-Line Support** – Supports options like `--debug` and `--file`.  
♦️- **Git** – For version control → [Git Documentation](https://git-scm.com/doc)

---

## 🛠️ **Technologies stack**
- **Python 3.9+** – Programming language  
- **Biopython** – For connecting with the PubMed API  
- **Pandas** – For data handling and CSV generation  
- **Poetry** – For dependency management  
- **argparse** – For CLI handling  

---

## 📥**Installation and Setup Process:**
### **1. Clone the Repository**
Open your terminal and run:
```bash
git clone https://github.com/your-username/research-paper-by-using-pubmed-api.git
cd research-paper-by-using-pubmed-api
```

### **2. Install Poetry**
Open your terminal and run:
```bash
pip install poetry
```
### **3. Initialize the Project with Poetry**
Set up the environment using Poetry:
```bash
poetry install
```

### **4. Set Up Python Environment (if needed)**
If you have multiple Python versions installed:
```bash
poetry env use python3.10
```
---

##🧑‍🔧**Configuration:**
Make sure to set your email in the script to comply with PubMed API guidelines:
```bash
Entrez.email = "your-email@example.com"
```
---

## **💡How to Use:**
### **1. Basic Search**
```bash
poetry run get-papers-list "cancer"
```
---

### **2. Save Output to a CSV File**
To save the output to a CSV file:
```bash
poetry run get-papers-list "AI in medicine" -m 20 -f results.csv
```

### **3. Enable Debug Mode**
To enable debug mode and print detailed output:
```bash
poetry run get-papers-list "drug discovery" -d
```
---

## **🔗 Project Structure**
```bash
Edit
├── fetch_pubmed.py        # Main script for fetching and processing papers
├── pyproject.toml         # Poetry project configuration file
├── README.md              # Project documentation
├── .gitignore             # Ignore unnecessary files
└── results.csv            # Output file (if specified)
```

----

## 📝 **Example CSV Output**
| Pubmed ID | Title | Publication Date | Non-academic Authors | Company Affiliations | Corresponding Author Email |
|-----------|-------|------------------|-----------------------|----------------------|----------------------------|
| 12345678  | AI in Drug Discovery | 2025 | John Doe, Jane Smith | XYZ Pharma, ABC Biotech | mkaif0262@gmail.com |

----

##🔧**Troubleshooting:**
### **1. Permission Denied for Poetry Commands
Try running commands with sudo:
```bash
sudo poetry install
```

### **2. Encoding Issues:**
If you get encoding errors, convert files to UTF-8:
```bash
iconv -f WINDOWS-1252 -t UTF-8 README.md > README-converted.md
```







