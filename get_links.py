# get_links.py
import urllib.request
import re

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=5431503&retmode=xml'
try:
    print("Fetching PMC XML...")
    xml = urllib.request.urlopen(url).read().decode('utf-8')
    print("Parsing links...")
    links = re.findall(r'href="([^"]+)"', xml)
    for link in links:
        print(link)
    
    # Also search for <supplementary-material> tags
    materials = re.findall(r'<supplementary-material[^>]*>(.*?)</supplementary-material>', xml, re.DOTALL)
    print("\nSupplementary Materials tags:")
    for mat in materials:
        print(mat.strip())
except Exception as e:
    print(f"Error: {e}")
