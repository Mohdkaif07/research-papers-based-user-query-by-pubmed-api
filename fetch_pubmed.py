from Bio import Entrez
import pandas as pd
import argparse

Entrez.email = "mkaif0262@gmail.com"

# Keywords for biotech and non-academic affiliation checks
BIOTECH_KEYWORDS = [
    "pharma", "biotech", "biopharma", "pharmaceutical", "biosciences",
    "med", "medical", "health", "therapeutics", "genomics", "life sciences",
    "drug", "biotechnology", "clinical", "biopharmaceutical"
]

NON_ACADEMIC_KEYWORDS = [
    "pharma", "biotech", "biopharma", "pharmaceutical", "biosciences",
    "med", "medical", "health", "therapeutics", "genomics", "life sciences",
    "drug", "biotechnology", "clinical", "biopharmaceutical",
    "labs", "corporation", "company", "LLC", "inc", "hospital"
]

# Check for biotech affiliation
def is_biotech(affiliation):
    return any(keyword in affiliation.lower() for keyword in BIOTECH_KEYWORDS) if affiliation else False

#  Check for non-academic affiliation
def is_non_academic(affiliation):
    return any(keyword in affiliation.lower() for keyword in NON_ACADEMIC_KEYWORDS) if affiliation else False

# Fetch papers from PubMed
def fetch_papers(query, max_results=20, debug=False):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        paper_ids = record.get("IdList", [])
        if not paper_ids:
            print("No papers found for the given query.")
            return []

        if debug:
            print(f"if not Found {len(paper_ids)} papers")

        #  Fetch papers in medline format
        handle = Entrez.efetch(db="pubmed", id=",".join(paper_ids), rettype="medline", retmode="text")
        records = handle.read()
        handle.close()

        return parse_medline_records(records, debug)

    except Exception as e:
        print(f"❌ Error during fetching: {e}")
        return []

#  Parse Medline records
def parse_medline_records(medline_records, debug):
    papers = []
    records = medline_records.strip().split("\n\n")

    for record in records:
        paper_data = {
            "Title": None,
            "Authors": [],
            "Affiliations": [],
            "Corresponding Author Email": None,
            "Pubmed ID": None,
            "Publication Date": None
        }

        lines = record.split("\n")
        for line in lines:
            if line.startswith("PMID- "):
                paper_data["Pubmed ID"] = line[6:].strip()
            elif line.startswith("TI  - "):
                paper_data["Title"] = line[6:].strip()
            elif line.startswith("DP  - "):
                paper_data["Publication Date"] = line[6:].strip()
            elif line.startswith("AU  - "):
                author = line[6:].strip()
                paper_data["Authors"].append(author)
            elif line.startswith("AD  - "):
                affiliation = line[6:].strip()
                paper_data["Affiliations"].append(affiliation)
            elif line.startswith("CA  - "):
                paper_data["Corresponding Author Email"] = line[6:].strip()

        if should_include_paper(paper_data, debug):
            papers.append(paper_data)

    return papers

#  Check if paper should be included based on biotech and non-academic rules
def should_include_paper(paper, debug):
    matched = False
    all_non_academic = True

    for affiliation in paper["Affiliations"]:
        if is_biotech(affiliation):
            matched = True

        if not is_non_academic(affiliation):
            all_non_academic = False

    if matched and all_non_academic:
        if debug:
            print(f"✅ Paper included: {paper['Title']}")
        return True
    else:
        if debug:
            print(f"❌ Paper excluded: {paper['Title']}")
        return False
    
#  Save to CSV
def convert_to_csv(papers, filename):
    if papers:
        df = pd.DataFrame(papers)
        df.to_csv(filename, index=False)
        print(f"✅ Results saved to {filename}")
    else:
        print("❌ No matching papers found.")

# Print to console
def print_to_console(papers):
    if papers:
        df = pd.DataFrame(papers)
        # Clean up the output by flattening lists for authors and affiliations
        df['Authors'] = df['Authors'].apply(lambda authors: ', '.join(authors))
        df['Affiliations'] = df['Affiliations'].apply(lambda affils: ', '.join(affils))
        df['Corresponding Author Email'] = df['Corresponding Author Email'].fillna("Not provided")
        
        # This will print the DataFrame to the console
        print(df.to_string(index=False))
    else:
        print("❌ No matching papers found.")

# Main function (CLI Handling)
def main():
    parser = argparse.ArgumentParser(description="Fetch research papers from PubMed")

    parser.add_argument("query", type=str, help="Search query for PubMed")
    parser.add_argument("-m", "--max_results", type=int, default=50, help="Maximum number of results to fetch")
    parser.add_argument("-f", "--file", type=str, help="File to save results (CSV). If not provided, output will be printed to console")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")
    args = parser.parse_args()

    papers = fetch_papers(args.query, max_results=args.max_results, debug=args.debug)

    if papers:
        if args.file:
            convert_to_csv(papers, args.file)
        else:
            print_to_console(papers)  # This will print the results to the terminal
    else:
        print("No matching papers found.")

if __name__ == "__main__":
    main()
