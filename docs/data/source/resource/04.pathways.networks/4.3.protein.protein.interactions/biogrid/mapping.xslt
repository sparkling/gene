<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  BioGRID -> Unified PPI Schema Mapping
  ============================================================
  Source: ./schema.json (BioGRID Tab 3.0 format)
  Target: ../schema.json (Unified 4.3 PPI Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ BioGRID_Interaction_ID      │ biogrid_interaction_id  │ Direct       │
  │ BioGRID_ID_A                │ biogrid_id_a            │ Direct       │
  │ BioGRID_ID_B                │ biogrid_id_b            │ Direct       │
  │ Entrez_Gene_A               │ interactor_a_id         │ Direct       │
  │ Entrez_Gene_B               │ interactor_b_id         │ Direct       │
  │ Official_Symbol_A           │ gene_symbol_a           │ Direct       │
  │ Official_Symbol_B           │ gene_symbol_b           │ Direct       │
  │ Systematic_Name_A           │ systematic_name_a       │ Direct       │
  │ Synonyms_A                  │ synonyms_a              │ Split array  │
  │ Experimental_System         │ experimental_system     │ Direct       │
  │ Experimental_System_Type    │ experimental_system_type│ Direct       │
  │ Throughput                  │ throughput              │ Direct       │
  │ Modification                │ modification            │ Direct       │
  │ Ontology_Term_IDs           │ ontology_term_ids       │ Split array  │
  │ SWISS-PROT_Accession_A      │ swissprot_accession_a   │ Direct       │
  │ Organism_ID_A               │ organism_a              │ Direct       │
  │ Organism_ID_B               │ organism_b              │ Direct       │
  │ Pubmed_ID                   │ publication_ids         │ Split array  │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - "-" -> null
  - Empty strings "" -> null

  NOTES:
  - BioGRID uses pipe (|) delimiter for multi-value fields
  - Experimental_System_Type is "physical" or "genetic"
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== REQUIRED FIELDS ========== -->

      <!-- Interactor A ID (Entrez Gene) -->
      <fn:string key="interactor_a_id">
        <xsl:value-of select="fn:string[@key='Entrez_Gene_A']"/>
      </fn:string>

      <!-- Interactor B ID (Entrez Gene) -->
      <fn:string key="interactor_b_id">
        <xsl:value-of select="fn:string[@key='Entrez_Gene_B']"/>
      </fn:string>

      <!-- ========== GENE SYMBOLS ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='Official_Symbol_A'] and fn:string[@key='Official_Symbol_A'] != '' and fn:string[@key='Official_Symbol_A'] != '-'">
          <fn:string key="gene_symbol_a">
            <xsl:value-of select="fn:string[@key='Official_Symbol_A']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbol_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:string[@key='Official_Symbol_B'] and fn:string[@key='Official_Symbol_B'] != '' and fn:string[@key='Official_Symbol_B'] != '-'">
          <fn:string key="gene_symbol_b">
            <xsl:value-of select="fn:string[@key='Official_Symbol_B']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_symbol_b"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BIOGRID-SPECIFIC IDS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='BioGRID_Interaction_ID']">
          <fn:number key="biogrid_interaction_id">
            <xsl:value-of select="fn:number[@key='BioGRID_Interaction_ID']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="biogrid_interaction_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='BioGRID_ID_A']">
          <fn:number key="biogrid_id_a">
            <xsl:value-of select="fn:number[@key='BioGRID_ID_A']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="biogrid_id_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='BioGRID_ID_B']">
          <fn:number key="biogrid_id_b">
            <xsl:value-of select="fn:number[@key='BioGRID_ID_B']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="biogrid_id_b"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SYSTEMATIC NAMES AND SYNONYMS ========== -->

      <xsl:choose>
        <xsl:when test="fn:string[@key='Systematic_Name_A'] and fn:string[@key='Systematic_Name_A'] != '' and fn:string[@key='Systematic_Name_A'] != '-'">
          <fn:string key="systematic_name_a">
            <xsl:value-of select="fn:string[@key='Systematic_Name_A']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="systematic_name_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Synonyms (pipe-delimited) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Synonyms_A'] and fn:string[@key='Synonyms_A'] != '' and fn:string[@key='Synonyms_A'] != '-'">
          <fn:array key="synonyms_a">
            <xsl:for-each select="tokenize(fn:string[@key='Synonyms_A'], '\|')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="synonyms_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EXPERIMENTAL DETAILS ========== -->

      <!-- Experimental system (e.g., Two-hybrid, Affinity Capture-MS) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Experimental_System'] and fn:string[@key='Experimental_System'] != ''">
          <fn:string key="experimental_system">
            <xsl:value-of select="fn:string[@key='Experimental_System']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="experimental_system"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Experimental system type (physical/genetic) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Experimental_System_Type'] and fn:string[@key='Experimental_System_Type'] != ''">
          <fn:string key="experimental_system_type">
            <xsl:value-of select="fn:string[@key='Experimental_System_Type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="experimental_system_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Throughput (high/low) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Throughput'] and fn:string[@key='Throughput'] != ''">
          <fn:string key="throughput">
            <xsl:value-of select="fn:string[@key='Throughput']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="throughput"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Modification -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Modification'] and fn:string[@key='Modification'] != '' and fn:string[@key='Modification'] != '-'">
          <fn:string key="modification">
            <xsl:value-of select="fn:string[@key='Modification']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="modification"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ONTOLOGY AND CROSS-REFS ========== -->

      <!-- Ontology terms (pipe-delimited) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Ontology_Term_IDs'] and fn:string[@key='Ontology_Term_IDs'] != '' and fn:string[@key='Ontology_Term_IDs'] != '-'">
          <fn:array key="ontology_term_ids">
            <xsl:for-each select="tokenize(fn:string[@key='Ontology_Term_IDs'], '\|')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ontology_term_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Swiss-Prot accession -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='SWISS-PROT_Accession_A'] and fn:string[@key='SWISS-PROT_Accession_A'] != '' and fn:string[@key='SWISS-PROT_Accession_A'] != '-'">
          <fn:string key="swissprot_accession_a">
            <xsl:value-of select="fn:string[@key='SWISS-PROT_Accession_A']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="swissprot_accession_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ORGANISM INFO ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='Organism_ID_A']">
          <fn:number key="organism_a">
            <xsl:value-of select="fn:number[@key='Organism_ID_A']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='Organism_ID_B']">
          <fn:number key="organism_b">
            <xsl:value-of select="fn:number[@key='Organism_ID_B']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_b"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PUBLICATIONS ========== -->

      <!-- PubMed IDs (pipe-delimited) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='Pubmed_ID'] and fn:string[@key='Pubmed_ID'] != '' and fn:string[@key='Pubmed_ID'] != '-'">
          <fn:array key="publication_ids">
            <xsl:for-each select="tokenize(fn:string[@key='Pubmed_ID'], '\|')">
              <fn:string><xsl:value-of select="normalize-space(.)"/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="publication_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">BioGRID</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='BioGRID_Interaction_ID']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
