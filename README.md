# 🧬 GENOMIC POWER - Hemoglobin Analysis Pipeline

A comprehensive genomic analysis pipeline for hemoglobin variants, focusing on pathogenicity prediction and clinical interpretation using machine learning and real clinical data.

## 📋 Overview

This project implements a complete workflow for analyzing hemoglobin gene variants, from raw ClinVar data processing to machine learning-based pathogenicity prediction. The pipeline integrates multiple genomic databases and APIs to provide robust variant classification.

## 🏗️ Architecture

### Core Components

1. **Data Processing Pipeline**
   - ClinVar VCF processing and annotation extraction
   - VEP (Variant Effect Predictor) data integration
   - Training set creation with labeled pathogenic/benign variants

2. **Genomic Mapping Engine**
   - Robust VEP CSQ field parsing with multiple format detection
   - MANE transcript prioritization
   - gnomAD integration with caching and exponential backoff
   - Reference genome handling (GRCh38/primary_assembly)

3. **Machine Learning Pipeline**
   - Feature engineering from genomic and clinical data
   - Pathogenicity prediction models
   - Cross-validation and performance evaluation

4. **Clinical Integration**
   - Mutation analysis and effect prediction
   - Clinical database integration
   - Report generation and visualization

## 📊 Data Sources

- **ClinVar**: Real clinical variant annotations
- **Ensembl**: Transcript and protein mapping
- **gnomAD**: Population frequency data
- **VEP**: Variant effect predictions
- **GenBank**: Gene sequences and annotations

## 🛠️ Technologies Used

### Programming Languages
- **Python 3.8+**: Core pipeline implementation

### Key Libraries
- **cyvcf2**: High-performance VCF parsing
- **pandas**: Data manipulation and analysis
- **scikit-learn**: Machine learning algorithms
- **pysam**: FASTA/sequence handling
- **requests**: API integrations
- **BioPython**: Biological sequence analysis

### System Tools
- **bcftools**: VCF manipulation
- **htslib**: HTS file handling
- **Mermaid**: Architecture visualization

## 🚀 Quick Start

### Prerequisites
```bash
# Install Python dependencies
pip install cyvcf2 pandas scikit-learn pysam requests biopython

# Install system tools (WSL/Linux)
sudo apt-get install bcftools htslib
```

### Basic Usage
```python
from hemoglobin_analysis.scripts.robust_genomic_mapper import RobustGenomicMapper

# Initialize mapper
mapper = RobustGenomicMapper(vep_data_path="data/processed/vep_extracted.tsv")

# Query variant information
result = mapper.get_variant_summary("HBB", "p.(Glu6Val)")
print(result)
```

## 📁 Project Structure

```
hemoglobin_analysis/
├── data/
│   ├── raw/                 # Raw ClinVar and reference data
│   ├── processed/           # Processed datasets and features
│   ├── hemoglobin_genes/    # Gene-specific data and sequences
│   └── cache/              # API response caching
├── scripts/
│   ├── build_training_set.py        # ClinVar processing
│   ├── robust_genomic_mapper.py     # Core mapping engine
│   ├── mutational_analysis.py       # Clinical integration
│   └── test_gnomad_cache.py         # Testing utilities
├── mutation_analysis/      # Analysis results and reports
└── docs/                   # Documentation and schemas
```

## 🔬 Key Features

### Robust Genomic Mapping
- **Multi-format CSQ parsing**: Handles various VEP output formats
- **MANE prioritization**: Selects clinically relevant transcripts
- **gnomAD integration**: Population frequency with smart caching
- **Reference handling**: Supports GRCh38 primary assembly

### Clinical Data Processing
- **Real ClinVar data**: Uses actual clinical annotations
- **Pathogenicity labels**: Automated benign/pathogenic classification
- **Feature engineering**: Comprehensive variant descriptors

### Performance Optimizations
- **Caching system**: Reduces API calls and improves speed
- **Batch processing**: Efficient handling of large datasets
- **Error recovery**: Robust handling of API failures

## 📈 Results

The pipeline successfully processes hemoglobin variants with:
- **High accuracy**: Clinically validated predictions
- **Comprehensive coverage**: Multiple genes and variant types
- **Scalable architecture**: Handles large genomic datasets
- **Production ready**: Robust error handling and logging

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **ClinVar** for clinical variant annotations
- **gnomAD** for population frequency data
- **Ensembl** for genomic annotations
- **VEP** for variant effect predictions

## 📞 Contact

**Project Author**: [agdiegorusso-stack](https://github.com/agdiegorusso-stack)

---

*Built with ❤️ for advancing genomic medicine and hemoglobinopathy research*</content>
<parameter name="filePath">c:\Users\Miste\Downloads\GENOMIC POWER\README.md
