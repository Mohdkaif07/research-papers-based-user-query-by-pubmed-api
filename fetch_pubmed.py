from Bio import Entrez
import pandas as pd
import time
import argparse

Entrez.email = "mkaif0262@gmail.com"

def is_biotech(affiliation):
    keywords =  [
        "pharma", "biotech", "biopharma", "pharmaceutical", "biosciences",
        "med", "medical", "health", "therapeutics", "genomics", "life sciences",
        "drug", "biotechnology", "clinical", "biopharmaceutical"
    ]
    if affiliation:
        return any(keyword in affiliation.lower() for keyword in keywords)
    return False    


def fetch_papers(query, max_result=10, debug=False):
    handle = Entrez.esearch (db = "pubmed", term = query, retmax = max_result)
    record = Entrez.read(handle)
    handle.close()

    paper_ids = record["IdList"]
    if debug:
        print(f"Found {len(paper_ids)} papers")

        
    papers = []

    for paper_id in paper_ids:
        handle = Entrez.efetch(db="pubmed", id=paper_id, rettype= "xml", retmode="text")
        paper = Entrez.read(handle)
        handle.close()
        
        article = paper["PubmedArticle"][0]["MedlineCitation"]["Article"]
        title = article.get("ArticleTitle", "N/A")
        year = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "N/A")
        authors = article.get("AuthorList", [])

        matched = False
        author_list = []
        non_academic_authors = []
        company_affiliations = []
        corresponding_author_email = "N/A"

        
        for author in authors:
            if "LastName" in author and "ForeName" in author:
                name = f"{author['ForeName']} {author['LastName']}"
                author_list.append(name)
                
                affiliation = ""
                if "AffiliationInfo" in author and author["AffiliationInfo"]:
                    affiliation = author["AffiliationInfo"][0].get("Affiliation", "")
                                     
                    if is_biotech(affiliation):
                        matched = True
                        non_academic_authors.append(name)
                        company_affiliations.append(affiliation)

                if "Identifier" in author:
                    email = author.get("Identifier", "")
                    if "@" in email:
                        corresponding_author_email = email

                if debug:
                    print(f"Author: {name}, Affiliation: {affiliation}")

                    
                        
        if matched:
         #     papers.append({
         #     "Title":title,
         #     "Authors": ",".join(author_list),
         #     "Year" : year,
         #     "Pubmed ID" : paper_id,
         #     "URL": f"https://pubmed.ncbi.nlm.nih.gov/{paper_id}/"
         # })
           papers.append({
                "Pubmed ID": paper_id,
                "Title": title,
                "Publication Date": year,
                "Non-academic Authors": ", ".join(non_academic_authors) if non_academic_authors else "N/A",
                "Company Affiliations": ", ".join(company_affiliations) if company_affiliations else "N/A",
                "Corresponding Author Email": corresponding_author_email,
                # "URL": f"https://pubmed.ncbi.nlm.nih.gov/{paper_id}/"
            })

        time.sleep(0.5)

    return papers

def convert_to_csv(papers, filename):
    if papers:
        df = pd.DataFrame(papers)
        df.to_csv(filename, index=False)
        print(f"Result saved to {filename}")
    else:
        print("No papers found with pharma")

# query = input("Enter search query: ")
# papers = fetch_papers("query", max_result=20)
# convert_to_csv(papers, 'pubmed_research_papers.csv')

def main():
    parser = argparse.ArgumentParser(description="Fetch research papers from PubMed")
    
    parser.add_argument("query", type=str, help="Search query for PubMed")
    parser.add_argument("-m", "--max_results", type=int, default=10, help="Maximum number of results to fetch")
    parser.add_argument("-f", "--file", type=str, help="File to save results (CSV)")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")

    args = parser.parse_args()

    # Fetch papers using the provided query
    papers = fetch_papers(args.query, max_result=args.max_results, debug=args.debug)

    if args.file:
        convert_to_csv(papers, args.file)
    else:
        for paper in papers:
            print(paper)

if __name__ == "__main__":
    main()    

# query = "cancer"
# max_results = 20
# file = "pubmed_results.csv"
# debug = True

# # Fetch and save papers
# papers = fetch_papers(query, max_result=max_results, debug=debug)
# convert_to_csv(papers, file)
        