#!/usr/bin/env bash
{
	test blastn && echo “blastn correctly installed”
} || {
	echo “blastn either not installed or not on the .bashrc”
}

{
	test mafft && echo “mafft correctly installed”
} || {
	echo “mafft either not installed or not on your bash profile”
}

{ #try direct
	test edirect && echo “edirect correctly installed”
} || {
	echo “entrez utilities not correctly installed”
}

{ #try tbl2asn
	test tbl2asn && echo “tbl2asn correctly installed”
} || {
	echo “tbl2asn not correctly installed”
}