# remove preamble

file=intro_bnp_density.Rmd

sed -i "s/^end.rcode.*/\`\`\`/" ${file}
sed -i "s/^<\!--begin.*/\`\`\`{r}/" ${file}

sed -i "s/^<\/body>//" ${file}
sed -i "s/^<\/html>//" ${file}

sed -i "s/^\[latex\]/\$\$/" ${file}
sed -i "s/^\[\/latex\]/\$\$/" ${file}

sed -i "s/\[latex\]\s*/\$/g" ${file}
sed -i "s/\s*\[\/latex\]/\$/g" ${file}
