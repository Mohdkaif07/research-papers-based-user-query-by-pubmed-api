# 🚀 Research Paper Fetching using PubMed API

This project is a Python-based CLI tool that allows users to search for research papers from the **PubMed database** based on a user-specified query. It filters the results to include only papers with at least one author affiliated with a **pharmaceutical** or **biotech** company and saves the output to a CSV file.

---

## 📌 **Key Features**
✅ **User-Specified Query** – Search for papers using any keyword or topic.  
✅ **PubMed API Integration** – Fetch research papers directly from PubMed using `Entrez`.  
✅ **Author Affiliation Filtering** – Filters papers based on pharma/biotech-related affiliations.  
✅ **Structured CSV Output** – Saves results to a CSV file with key details.  
✅ **Command-Line Support** – Supports options like `--debug` and `--file`.  

---

## 🛠️ **Technologies Used**
- **Python 3.9+** – Programming language  
- **Biopython** – For connecting with the PubMed API  
- **Pandas** – For data handling and CSV generation  
- **Poetry** – For dependency management  
- **argparse** – For CLI handling  

---

## 🚀 **Installation and Setup**
### **1. Clone the Repository**
Open your terminal and run:
```bash
git clone https://github.com/your-username/research-paper-by-using-pubmed-api.git
cd research-paper-by-using-pubmed-api
