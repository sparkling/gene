<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  SWISS-MODEL -> Unified Protein Structures Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Structures Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field            │ Target Field        │ Transform    │
  ├─────────────────────────┼─────────────────────┼──────────────┤
  │ model_id                │ entry_id            │ Direct       │
  │ uniprot_ac              │ uniprot_accession   │ Direct       │
  │ from/to                 │ model_range         │ Restructure  │
  │ template                │ template            │ Direct       │
  │ identity                │ template_identity   │ Direct       │
  │ coverage                │ coverage            │ Direct       │
  │ qmean                   │ qmean               │ Direct       │
  │ qmean_disco             │ qmean_disco         │ Direct       │
  │ (computed)              │ structure_type      │ Constant     │
  │ (computed)              │ _source             │ Metadata     │
  └─────────────────────────┴─────────────────────┴──────────────┘

  NOTES:
  - SWISS-MODEL provides homology models (predicted_homology type)
  - Models flattened from nested array structure
  - Quality metrics mapped to unified schema fields
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    exclude-result-prefixes="xs fn map array">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- Process each model within a SWISS-MODEL entry -->
  <xsl:template match="fn:map">
    <xsl:variable name="uniprot" select="fn:string[@key='uniprot_ac']"/>

    <fn:array>
      <xsl:for-each select="fn:array[@key='models']/fn:map">
        <fn:map>
          <!-- Entry ID from model_id -->
          <fn:string key="entry_id">
            <xsl:value-of select="fn:string[@key='model_id']"/>
          </fn:string>

          <!-- UniProt accession -->
          <fn:string key="uniprot_accession">
            <xsl:value-of select="$uniprot"/>
          </fn:string>

          <!-- Structure type: homology model -->
          <fn:string key="structure_type">predicted_homology</fn:string>

          <!-- Template information -->
          <fn:string key="template">
            <xsl:value-of select="fn:string[@key='template']"/>
          </fn:string>

          <!-- Template identity -->
          <fn:number key="template_identity">
            <xsl:value-of select="fn:number[@key='identity']"/>
          </fn:number>

          <!-- Template similarity -->
          <xsl:if test="fn:number[@key='similarity']">
            <fn:number key="template_similarity">
              <xsl:value-of select="fn:number[@key='similarity']"/>
            </fn:number>
          </xsl:if>

          <!-- Coverage -->
          <fn:number key="coverage">
            <xsl:value-of select="fn:number[@key='coverage']"/>
          </fn:number>

          <!-- Model range -->
          <fn:map key="model_range">
            <fn:number key="from">
              <xsl:value-of select="fn:number[@key='from']"/>
            </fn:number>
            <fn:number key="to">
              <xsl:value-of select="fn:number[@key='to']"/>
            </fn:number>
          </fn:map>

          <!-- QMEAN score -->
          <fn:number key="qmean">
            <xsl:value-of select="fn:number[@key='qmean']"/>
          </fn:number>

          <!-- QMEANDisCo score -->
          <xsl:if test="fn:number[@key='qmean_disco']">
            <fn:number key="qmean_disco">
              <xsl:value-of select="fn:number[@key='qmean_disco']"/>
            </fn:number>
          </xsl:if>

          <!-- Oligomeric state -->
          <xsl:if test="fn:string[@key='oligo_state']">
            <fn:string key="oligo_state">
              <xsl:value-of select="fn:string[@key='oligo_state']"/>
            </fn:string>
          </xsl:if>

          <!-- Model created date -->
          <xsl:if test="fn:string[@key='created_date']">
            <fn:string key="model_created_date">
              <xsl:value-of select="fn:string[@key='created_date']"/>
            </fn:string>
          </xsl:if>

          <!-- Source metadata -->
          <fn:map key="_source">
            <fn:string key="database">SWISS-MODEL</fn:string>
            <fn:string key="original_id">
              <xsl:value-of select="fn:string[@key='model_id']"/>
            </fn:string>
            <fn:string key="uniprot_ac">
              <xsl:value-of select="$uniprot"/>
            </fn:string>
          </fn:map>

        </fn:map>
      </xsl:for-each>
    </fn:array>
  </xsl:template>

</xsl:stylesheet>
