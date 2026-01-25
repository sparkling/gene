<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  RefSeq -> Unified Protein Sequences Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Protein Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ accession               │ accession       │ Direct       │
  │ definition              │ definition      │ Direct       │
  │ organism                │ organism        │ Direct       │
  │ taxonomy                │ taxonomy        │ Direct       │
  │ gene                    │ gene            │ Direct       │
  │ gene_id                 │ gene_id         │ Direct       │
  │ sequence                │ sequence        │ Direct       │
  │ length                  │ length          │ Direct       │
  │ features                │ features        │ Restructure  │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - RefSeq is a core NCBI resource
  - Direct mapping to unified protein schema
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

  <xsl:template match="fn:map">
    <fn:map>
      <!-- Primary accession -->
      <fn:string key="accession">
        <xsl:value-of select="fn:string[@key='accession']"/>
      </fn:string>

      <!-- Sequence -->
      <fn:string key="sequence">
        <xsl:value-of select="fn:string[@key='sequence']"/>
      </fn:string>

      <!-- Length -->
      <fn:number key="length">
        <xsl:value-of select="fn:number[@key='length']"/>
      </fn:number>

      <!-- Organism -->
      <fn:string key="organism">
        <xsl:value-of select="fn:string[@key='organism']"/>
      </fn:string>

      <!-- Taxonomy -->
      <xsl:if test="fn:array[@key='taxonomy']">
        <fn:array key="taxonomy">
          <xsl:for-each select="fn:array[@key='taxonomy']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Gene -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='gene']">
          <fn:string key="gene">
            <xsl:value-of select="fn:string[@key='gene']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Gene ID -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='gene_id']">
          <fn:number key="gene_id">
            <xsl:value-of select="fn:number[@key='gene_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Definition -->
      <fn:string key="definition">
        <xsl:value-of select="fn:string[@key='definition']"/>
      </fn:string>

      <!-- Features -->
      <xsl:if test="fn:array[@key='features']">
        <fn:array key="features">
          <xsl:for-each select="fn:array[@key='features']/fn:map">
            <fn:map>
              <fn:string key="type">
                <xsl:value-of select="fn:string[@key='key']"/>
              </fn:string>
              <fn:map key="location">
                <xsl:variable name="loc" select="fn:string[@key='location']"/>
                <xsl:analyze-string select="$loc" regex="^([0-9]+)\.\.([0-9]+)$">
                  <xsl:matching-substring>
                    <fn:number key="start"><xsl:value-of select="regex-group(1)"/></fn:number>
                    <fn:number key="end"><xsl:value-of select="regex-group(2)"/></fn:number>
                  </xsl:matching-substring>
                </xsl:analyze-string>
              </fn:map>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Cross-references -->
      <xsl:if test="fn:array[@key='dbxrefs']">
        <fn:array key="dbxrefs">
          <xsl:for-each select="fn:array[@key='dbxrefs']/fn:map">
            <fn:map>
              <fn:string key="database">
                <xsl:value-of select="fn:string[@key='database']"/>
              </fn:string>
              <fn:string key="id">
                <xsl:value-of select="fn:string[@key='id']"/>
              </fn:string>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- RefSeq-specific fields -->
      <xsl:if test="fn:string[@key='coded_by']">
        <fn:string key="coded_by">
          <xsl:value-of select="fn:string[@key='coded_by']"/>
        </fn:string>
      </xsl:if>

      <xsl:if test="fn:string[@key='chromosome']">
        <fn:string key="chromosome">
          <xsl:value-of select="fn:string[@key='chromosome']"/>
        </fn:string>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">RefSeq</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='accession']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
