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
    if affiliation:
        # Explicitly check if the affiliation includes known academic institutions or university-related terms
        academic_keywords = ["university", "school of", "college", "institute"]
        # If it's an academic affiliation, return False
        if any(keyword in affiliation.lower() for keyword in academic_keywords):
            return False
        return any(keyword in affiliation.lower() for keyword in NON_ACADEMIC_KEYWORDS)
    return False

# Fetch papers from PubMed
def fetch_papers(query, max_results=100, debug=False):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        paper_ids = record.get("IdList", [])
        if not paper_ids:
            print("No papers found for the given query.")
            return []

        if debug:
            print(f"Found {len(paper_ids)} papers")

        #  Fetch papers in medline format
        handle = Entrez.efetch(db="pubmed", id=",".join(paper_ids), rettype="medline", retmode="text")
        records = handle.read()
        handle.close()

        return parse_medline_records(records, debug)

    except Exception as e:
        print(f"❌ Error during fetching: {e}")
        return []

# Parse Medline records
def parse_medline_records(medline_records, debug):
    papers = []
    records = medline_records.strip().split("\n\n")

    for record in records:
        paper_data = {
            "Title": None,
            "Authors": [],
            "Affiliations": [],
            "Non-academic Author": [],
            "Company Affiliation": [],
            "Corresponding Author Email": "Not provided",  # Default to "Not provided"
            "Pubmed ID": None,
            "Publication Date": None  # Initialize with None
        }

        lines = record.split("\n")
        for line in lines:
            if line.startswith("PMID- "):
                paper_data["Pubmed ID"] = line[6:].strip()
            elif line.startswith("TI  - "):
                paper_data["Title"] = line[6:].strip()
            elif line.startswith("DP  - "):  # Publication Date
                paper_data["Publication Date"] = line[6:].strip()
            elif line.startswith("AU  - "):  # Authors
                author = line[6:].strip()
                paper_data["Authors"].append(author)
            elif line.startswith("AD  - "):  # Affiliations
                affiliation = line[6:].strip()
                paper_data["Affiliations"].append(affiliation)
                # Check if the author has a non-academic affiliation
                if is_non_academic(affiliation):
                    paper_data["Non-academic Author"].append(paper_data["Authors"][-1])  # Append last author
                # Check if affiliation is company-related
                if is_biotech(affiliation):
                    paper_data["Company Affiliation"].append(affiliation)
            elif line.startswith("CA  - "):  # Corresponding Author Email
                paper_data["Corresponding Author Email"] = line[6:].strip()

        # Handle case where no publication date is provided
        if not paper_data["Publication Date"]:
            paper_data["Publication Date"] = "Not provided"

        # If corresponding author email is not found, default to "Not provided"
        if not paper_data["Corresponding Author Email"]:
            paper_data["Corresponding Author Email"] = "Not provided"

        if should_include_paper(paper_data, debug):
            papers.append(paper_data)

    return papers

# Check if paper should be included based on biotech and non-academic rules
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

# Save to CSV
def convert_to_csv(papers, filename):
    if papers:
        # Reorder columns as required: Pubmed ID, Title, Publication Date, Non-academic Author, Company Affiliation, Corresponding Author Email
        ordered_papers = [{
            "Pubmed ID": paper["Pubmed ID"],
            "Title": paper["Title"],
            "Publication Date": paper["Publication Date"],
            "Non-academic Author": ', '.join(paper["Non-academic Author"]),
            "Company Affiliation": ', '.join(paper["Company Affiliation"]),
            "Corresponding Author Email": paper["Corresponding Author Email"]
        } for paper in papers]

        # Create DataFrame and save to CSV
        df = pd.DataFrame(ordered_papers)
        df.to_csv(filename, index=False)
        print(f"✅ Results saved to {filename}")
    else:
        print("❌ No matching papers found.")

# Print to console
def print_to_console(papers):
    if papers:
        # Reorder columns and clean up the output for console print
        ordered_papers = [{
            "Pubmed ID": paper["Pubmed ID"],
            "Title": paper["Title"],
            "Publication Date": paper["Publication Date"],
            "Non-academic Author": ', '.join(paper["Non-academic Author"]),
            "Company Affiliation": ', '.join(paper["Company Affiliation"]),
            "Corresponding Author Email": paper["Corresponding Author Email"]
        } for paper in papers]

        # Create DataFrame and print to console
        df = pd.DataFrame(ordered_papers)
        print(df.to_string(index=False))
    else:
        print("❌ No matching papers found.")

# Main function (CLI Handling)
def main():
    parser = argparse.ArgumentParser(description="Fetch research papers from PubMed")

    parser.add_argument("query", type=str, help="Search query for PubMed")
    parser.add_argument("-m", "--max_results", type=int, default=100, help="Maximum number of results to fetch")
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
