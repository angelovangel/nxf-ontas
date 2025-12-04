#!/usr/bin/env python3
"""
Generate an HTML report from coverage histogram files (.hist)
Usage: python bedtools-report.py sample1.hist [sample2.hist sample3.hist ...] -o output.html
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

def parse_hist_file(filepath):
    """Parse a .hist file and return list of records"""
    data = []
    with open(filepath, 'r') as f:
        # Skip the first line, assuming it is a header
        next(f, None)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                data.append({
                    'chr': parts[0],
                    'start': int(parts[1]) if parts[1].isdigit() else 0,
                    'end': int(parts[2]) if parts[2].isdigit() else 0,
                    'gene': parts[3],
                    'coverage': int(parts[4]),
                    'count': int(parts[5]),
                    'total': int(parts[6]),
                    'fraction': float(parts[7])
                })
    return data

def calculate_cumulative_coverage(gene_data):
    """Calculate cumulative coverage metrics"""
    total_bases = gene_data[0]['total'] if gene_data else 0
    
    # Calculate bases with coverage >= threshold
    bases_1x = sum(r['count'] for r in gene_data if r['coverage'] >= 1)
    bases_10x = sum(r['count'] for r in gene_data if r['coverage'] >= 10)
    bases_20x = sum(r['count'] for r in gene_data if r['coverage'] >= 20)
    bases_30x = sum(r['count'] for r in gene_data if r['coverage'] >= 30)
    
    return {
        'total': total_bases,
        'cov_1x': bases_1x,
        'cov_10x': bases_10x,
        'cov_20x': bases_20x,
        'cov_30x': bases_30x,
        'pct_1x': (bases_1x / total_bases * 100) if total_bases > 0 else 0,
        'pct_10x': (bases_10x / total_bases * 100) if total_bases > 0 else 0,
        'pct_20x': (bases_20x / total_bases * 100) if total_bases > 0 else 0,
        'pct_30x': (bases_30x / total_bases * 100) if total_bases > 0 else 0,
    }


def generate_html_report(samples_data, output_file):
    """Generate HTML report from multiple samples"""
    
    # Get unique genes
    all_genes = set()
    for sample_name, sample_data in samples_data.items():
        all_genes.update(r['gene'] for r in sample_data)
    genes = sorted(all_genes)
    
    # Pre-generate select options for the Regions filter (used by Select2)
    region_options = "".join([f'<option value="{gene.lower()}">{gene}</option>' for gene in genes])
    
    # Pre-generate select options for the Samples filter (used by Select2)
    sample_names = sorted(samples_data.keys())
    sample_options = "".join([f'<option value="{name.lower()}">{name}</option>' for name in sample_names])
    
    # Calculate total bases per region
    region_totals = {}
    for gene in genes:
        for sample_name, sample_data in samples_data.items():
            gene_data = [r for r in sample_data if r['gene'] == gene]
            if gene_data:
                region_totals[gene] = gene_data[0]['total']
                break
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Coverage Histogram Report</title>
  
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
  
  <link href="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/css/select2.min.css" rel="stylesheet" />
  
  <style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ 
      /* Use Inter as the main UI font (loaded via CDN) */
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
      background: #f5f5f5;
      min-height: 100vh;
      padding: 20px;
    }}
    .container {{
      max-width: 1800px;
      margin: 0 auto;
      background: white;
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.1);
      overflow: hidden;
    }}
    .header {{
      background: #60a5fa;
      color: white;
      padding: 20px;
      text-align: center;
    }}
    .header h1 {{ font-size: 1.5em; margin-bottom: 5px; }}
    .header p {{ opacity: 0.9; font-size: 0.9em; }}
    .content {{ padding: 20px; }}
    .stats {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 15px;
      margin-bottom: 20px;
    }}
    .stat-card {{
      background: #60a5fa;
      color: white;
      padding: 15px;
      border-radius: 8px;
      text-align: center;
    }}
    .stat-card h3 {{ font-size: 0.8em; opacity: 0.9; margin-bottom: 3px; }}
    .stat-card .value {{ font-size: 1.5em; font-weight: bold; }}
    .controls {{
      display: flex;
      gap: 20px;
      align-items: flex-end;
      margin-bottom: 20px;
    }}
    .filter-group {{
      flex-grow: 1;
    }}
    .filter-label {{
      display: block;
      margin-bottom: 5px;
      font-weight: 600;
      font-size: 0.9em;
    }}
    .filter-box {{
      padding: 10px;
      width: 100%;
      border: 2px solid #e2e8f0;
      border-radius: 6px;
      font-size: 16px;
    }}
    .filter-box:focus {{
      outline: none;
      border-color: #60a5fa;
    }}
    /* --- Select2 Custom Styling --- */
    .filter-group .select2-container {{
      width: 100% !important; 
    }}
    .select2-container .select2-selection--multiple {{
        min-height: 40px; 
        border: 2px solid #e2e8f0 !important;
        border-radius: 6px;
        padding: 5px 10px;
    }}
    .select2-container--default.select2-container--focus .select2-selection--multiple {{
        border-color: #60a5fa !important;
        outline: 0;
    }}
    .select2-container--default .select2-selection--multiple .select2-selection__choice {{
        background-color: #e0f2fe;
        border: 1px solid #90cdf4;
        color: #1e40af;
        padding: 3px 5px;
        border-radius: 4px;
    }}
    /* --- End Select2 Custom Styling --- */
    
    .table-container {{
      overflow-x: auto;
      border-radius: 8px;
      border: 1px solid #e2e8f0;
      max-height: 600px;
      overflow-y: auto;
      position: relative;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: white;
    }}
    thead {{
      position: sticky;
      top: 0;
      z-index: 100;
      background: #f8fafc;
    }}
    th, td {{
      padding: 10px;
      text-align: left;
      border-bottom: 1px solid #e2e8f0;
      font-size: 0.9em;
      /* Monospace font for data */
      font-family: 'Courier New', monospace; 
    }}
    th {{
      background: #f8fafc;
      font-weight: 600;
      /* Sans-serif for headers */
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
    }}
    tr:hover {{ background: #f8fafc; }}
    .coverage-cell {{
      text-align: right;
    }}
    .pct-bar {{
      display: inline-block;
      height: 16px;
      background: #60a5fa;
      border-radius: 2px;
      min-width: 2px;
      vertical-align: middle;
      margin-left: 5px;
    }}
    .sample-col {{
      font-weight: 600;
      background: #f8fafc;
    }}
    .export-btn {{
      background: #60a5fa;
      color: white;
      border: none;
      padding: 10px 20px;
      border-radius: 6px;
      cursor: pointer;
      font-size: 14px;
      font-weight: 600;
      white-space: nowrap; 
    }}
    .export-btn:hover {{
      background: #3b82f6;
    }}
    th.sortable {{
      cursor: pointer;
      user-select: none;
      position: relative;
      padding-right: 20px;
    }}
    th.sortable:hover {{
      background: #e2e8f0;
    }}
    th.sortable::after {{
      content: '⇅';
      position: absolute;
      right: 5px;
      opacity: 0.3;
      font-size: 0.8em;
    }}
    th.sortable.asc::after {{
      content: '▲';
      opacity: 1;
    }}
    th.sortable.desc::after {{
      content: '▼';
      opacity: 1;
    }}
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <h1>Coverage Histogram Report</h1>
      <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    <div class="content">
      <div class="stats">
        <div class="stat-card">
          <h3>Total Genes/Regions</h3>
          <div class="value">{len(genes)}</div>
        </div>
        <div class="stat-card">
          <h3>Samples</h3>
          <div class="value">{len(samples_data)}</div>
        </div>
        <div class="stat-card">
          <h3>Total Bases Analyzed</h3>
          <div class="value">{sum(region_totals.values()):,}</div>
        </div>
      </div>
      
      <div class="controls">
        <div class="filter-group">
          <label for="sampleFilter" class="filter-label">Filter Samples (Select/Search)</label>
          <select id="sampleFilter" multiple="multiple" style="width: 100%;">
            {sample_options}
          </select>
        </div>
        <div class="filter-group">
          <label for="regionFilter" class="filter-label">Filter Regions (Select/Search)</label>
          <select id="regionFilter" multiple="multiple" style="width: 100%;">
            {region_options}
          </select>
        </div>
        <button class="export-btn" onclick="exportToCSV()">Export CSV</button>
      </div>
      
      <div class="table-container">
        <table id="dataTable">
          <thead>
            <tr>
              <th rowspan="2" class="sortable" onclick="sortTable(0)">Sample</th>
              <th rowspan="2" class="sortable" onclick="sortTable(1)">Chr</th>
              <th rowspan="2" class="sortable" onclick="sortTable(2)">Gene/Region</th>
              <th rowspan="2" class="sortable" onclick="sortTable(3)">Location</th>
              <th rowspan="2" class="sortable" onclick="sortTable(4)">Total Bases</th>
              <th colspan="4" style="text-align: center; border-bottom: 1px solid #cbd5e1;">Cumulative Coverage</th>
            </tr>
            <tr>
              <th class="sortable" onclick="sortTable(5)" style="text-align: right;">≥1x (%)</th>
              <th class="sortable" onclick="sortTable(6)" style="text-align: right;">≥10x (%)</th>
              <th class="sortable" onclick="sortTable(7)" style="text-align: right;">≥20x (%)</th>
              <th class="sortable" onclick="sortTable(8)" style="text-align: right;">≥30x (%)</th>
            </tr>
          </thead>
          <tbody>
"""
    
    # Generate table rows
    for gene in genes:
        for sample_name, sample_data in samples_data.items():
            gene_data = [r for r in sample_data if r['gene'] == gene]
            
            if not gene_data:
                continue
                
            location = gene_data[0]
            cov_stats = calculate_cumulative_coverage(gene_data)
            
            # Row data attributes are used for sorting, filtering, and clean CSV export
            html += f"""
            <tr data-gene="{gene.lower()}" data-sample="{sample_name.lower()}" 
                data-chr="{location['chr']}"
                data-location="{location['start']:,}-{location['end']:,}"
                data-total="{cov_stats['total']}" 
                data-cov1="{cov_stats['pct_1x']:.1f}" 
                data-cov10="{cov_stats['pct_10x']:.1f}"
                data-cov20="{cov_stats['pct_20x']:.1f}"
                data-cov30="{cov_stats['pct_30x']:.1f}">
              <td class="sample-col">{sample_name}</td>
              <td>{location['chr']}</td>
              <td><strong>{gene}</strong></td>
              <td>{location['start']:,}-{location['end']:,}</td>
              <td style="text-align: right;">{cov_stats['total']:,}</td>
              <td class="coverage-cell">
                {cov_stats['pct_1x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_1x']*1.5}px;"></span>
              </td>
              <td class="coverage-cell">
                {cov_stats['pct_10x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_10x']*1.5}px;"></span>
              </td>
              <td class="coverage-cell">
                {cov_stats['pct_20x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_20x']*1.5}px;"></span>
              </td>
              <td class="coverage-cell">
                {cov_stats['pct_30x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_30x']*1.5}px;"></span>
              </td>
            </tr>
"""
    
    html += """
          </tbody>
        </table>
      </div>
    </div>
  </div>
  
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/js/select2.min.js"></script>
  
  <script>
    let sortDirection = {};
    
    // --- Select2 Initialization and Filter Binding ---
    $(document).ready(function() {
      // Initialize Select2 on the sample filter
      $('#sampleFilter').select2({
        placeholder: "Select or type sample names...",
        allowClear: true,
        closeOnSelect: false 
      }).on('change', filterTable); 
      
      // Initialize Select2 on the region filter
      $('#regionFilter').select2({
        placeholder: "Select or type regions...",
        allowClear: true,
        closeOnSelect: false 
      }).on('change', filterTable); 
      
      // Apply initial filter
      filterTable();
    });
    // -------------------------------------------------
    
    function exportToCSV() {
      const table = document.getElementById('dataTable');
      let csv = [];
      
      // Headers
      const headers = [
        'Sample', 
        'Chr', 
        'Gene/Region', 
        'Location', 
        'Total Bases', 
        'Pct_Ge_1x', 
        'Pct_Ge_10x', 
        'Pct_Ge_20x', 
        'Pct_Ge_30x'
      ];
      csv.push(headers.join(','));
      
      // Data rows (only visible ones)
      const dataRows = table.querySelectorAll('tbody tr');
      dataRows.forEach(row => {
        if (row.style.display !== 'none') {
          // Use data attributes for clean data extraction
          const cols = [
            row.getAttribute('data-sample'),
            row.getAttribute('data-chr'),
            row.getAttribute('data-gene'),
            row.getAttribute('data-location').replace(/,/g, ''), // remove commas from location
            row.getAttribute('data-total'),
            row.getAttribute('data-cov1'),
            row.getAttribute('data-cov10'),
            row.getAttribute('data-cov20'),
            row.getAttribute('data-cov30'),
          ];
          
          // Escape quotes and wrap in quotes if contains comma
          const safeCols = cols.map(text => {
            text = text.trim();
            if (text.includes(',') || text.includes('"')) {
              text = '"' + text.replace(/"/g, '""') + '"';
            }
            return text;
          });
          
          csv.push(safeCols.join(','));
        }
      });
      
      // Create download
      const csvContent = csv.join('\\n');
      const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
      const link = document.createElement('a');
      const url = URL.createObjectURL(blob);
      link.setAttribute('href', url);
      link.setAttribute('download', 'coverage_report_' + new Date().getTime() + '.csv');
      link.style.visibility = 'hidden';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }
    
    function sortTable(columnIndex) {
      const table = document.getElementById('dataTable');
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.querySelectorAll('tr'));
      
      // Toggle sort direction
      sortDirection[columnIndex] = sortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
      const isAsc = sortDirection[columnIndex] === 'asc';
      
      // Update header styling
      const headers = table.querySelectorAll('th.sortable');
      headers.forEach(h => {
        h.classList.remove('asc', 'desc');
      });
      headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');
      
      rows.sort((a, b) => {
        let aVal, bVal;
        
        if (columnIndex === 0) { // Sample
          aVal = a.getAttribute('data-sample');
          bVal = b.getAttribute('data-sample');
        } else if (columnIndex === 1) { // Chr
          aVal = a.getAttribute('data-chr');
          bVal = b.getAttribute('data-chr');
        } else if (columnIndex === 2) { // Gene
          aVal = a.getAttribute('data-gene');
          bVal = b.getAttribute('data-gene');
        } else if (columnIndex === 3) { // Location
          // Use the start coordinate for location sort
          const aLoc = a.getAttribute('data-location').split('-')[0].replace(/,/g, '');
          const bLoc = b.getAttribute('data-location').split('-')[0].replace(/,/g, '');
          aVal = parseInt(aLoc);
          bVal = parseInt(bLoc);
        } else if (columnIndex === 4) { // Total Bases
          aVal = parseFloat(a.getAttribute('data-total'));
          bVal = parseFloat(b.getAttribute('data-total'));
        } else if (columnIndex === 5) { // ≥1x
          aVal = parseFloat(a.getAttribute('data-cov1'));
          bVal = parseFloat(b.getAttribute('data-cov1'));
        } else if (columnIndex === 6) { // ≥10x
          aVal = parseFloat(a.getAttribute('data-cov10'));
          bVal = parseFloat(b.getAttribute('data-cov10'));
        } else if (columnIndex === 7) { // ≥20x
          aVal = parseFloat(a.getAttribute('data-cov20'));
          bVal = parseFloat(b.getAttribute('data-cov20'));
        } else if (columnIndex === 8) { // ≥30x
          aVal = parseFloat(a.getAttribute('data-cov30'));
          bVal = parseFloat(b.getAttribute('data-cov30'));
        }
        
        if (typeof aVal === 'string') {
          return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        } else {
          // Handle cases where aVal or bVal might be NaN/null for a safer comparison
          if (isNaN(aVal) || isNaN(bVal)) {
            return 0;
          }
          return isAsc ? aVal - bVal : bVal - aVal;
        }
      });
      
      rows.forEach(row => tbody.appendChild(row));
    }
    
    function filterTable() {
      // Get selected samples from Select2 (returns an array of values, or null if none selected)
      const selectedSamples = $('#sampleFilter').val();
      const sampleFilters = selectedSamples || [];
      
      // Get selected regions from Select2
      const selectedRegions = $('#regionFilter').val(); 
      const regionFilters = selectedRegions || [];
      
      const table = document.getElementById('dataTable');
      const rows = table.getElementsByTagName('tbody')[0].getElementsByTagName('tr');
      
      for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        const sample = row.getAttribute('data-sample');
        const gene = row.getAttribute('data-gene');
        
        // 1. Sample Match: If filters exist, the row sample must be in the selected list.
        let sampleMatch = true;
        if (sampleFilters.length > 0) {
            sampleMatch = sampleFilters.includes(sample);
        }
        
        // 2. Region Match: If filters exist, the row gene must be in the selected list.
        let regionMatch = true;
        if (regionFilters.length > 0) {
          regionMatch = regionFilters.includes(gene);
        }
        
        if (sampleMatch && regionMatch) {
          row.style.display = '';
        } else {
          row.style.display = 'none';
        }
      }
    }
  </script>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"Report generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Generate coverage histogram HTML report')
    parser.add_argument('samples', nargs='+', help='One or more .hist files')
    parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    
    args = parser.parse_args()
    
    samples_data = {}
    
    for sample_file in args.samples:
        sample_name = Path(sample_file).stem
        print(f"Loading {sample_file}...")
        samples_data[sample_name] = parse_hist_file(sample_file)
    
    print("Generating HTML report...")
    generate_html_report(samples_data, args.output)
    
    print(f"\nDone! Open {args.output} in your browser to view the report.")
    print(f"Total samples: {len(samples_data)}")
    print(f"Sample names: {', '.join(samples_data.keys())}")

if __name__ == "__main__":
    main()