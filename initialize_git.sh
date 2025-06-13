#!/bin/bash

# Initialize Git repository for Chimeric Detective
# Run this script after creating your GitHub repository

set -e

echo "🔧 Initializing Git repository for Chimeric Detective..."

# Initialize git repository if not already done
if [ ! -d ".git" ]; then
    echo "📁 Initializing Git repository..."
    git init
else
    echo "📁 Git repository already exists"
fi

# Add all files
echo "📝 Adding files to Git..."
git add .

# Create initial commit
echo "💾 Creating initial commit..."
git commit -m "Initial commit: Chimeric Detective v1.0.0

🔬 Comprehensive tool for detecting and resolving chimeric contigs in viral metagenomic assemblies

Features:
- Multi-method chimera detection (coverage, GC, k-mer, orientation)
- Intelligent classification (technical artifacts vs biological recombination)
- Automated resolution with precise breakpoint splitting
- Interactive HTML reports with visualizations
- Multi-sample processing with parallel execution
- Comprehensive CLI with flexible input support

Generated with Claude Code 🤖"

# Check if remote origin exists
if git remote get-url origin >/dev/null 2>&1; then
    echo "🔗 Remote origin already configured"
    git remote -v
else
    echo "❓ Remote origin not configured. Please set it up:"
    echo "   git remote add origin git@github.com:megjohnson1999/chimeric-contig-detector.git"
fi

echo ""
echo "✅ Git repository initialized successfully!"
echo ""
echo "🚀 Next steps:"
echo "1. Repository 'chimeric-contig-detector' should already exist on GitHub"
echo "2. Set the remote origin:"
echo "   git remote add origin git@github.com:megjohnson1999/chimeric-contig-detector.git"
echo "3. Push to GitHub:"
echo "   git branch -M main"
echo "   git push -u origin main"
echo ""
echo "📚 After pushing:"
echo "- Enable GitHub Actions for CI/CD"
echo "- Set up branch protection rules"
echo "- Configure Codecov for test coverage"
echo "- Add collaborators if needed"
echo ""
echo "🔬 Your Chimeric Detective repository is ready for GitHub!"