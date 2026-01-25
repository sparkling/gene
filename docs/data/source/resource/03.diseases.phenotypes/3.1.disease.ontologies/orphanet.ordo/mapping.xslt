<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Orphanet/ORDO -> Disease Ontologies Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Orphanet Rare Disease Ontology)
  Target: ../schema.json (3.1 Disease Ontologies Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ OrphaCode                   │ orpha_code              │ Direct       │
  │ Name                        │ name                    │ Direct       │
  │ Definition                  │ definition              │ Direct       │
  │ Synonym                     │ synonyms                │ Array        │
  │ OrphaGroup                  │ orpha_group             │ Enum         │
  │ DisorderType                │ disorder_type           │ Object       │
  │ AssociationType             │ association_type        │ Direct       │
  │ AssociationStatus           │ association_status      │ Enum         │
  │ ExternalReference           │ xref                    │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - Orphanet uses numeric OrphaCode
  - OrphaGroup: Disorder, Group, Subtype
  - Focus on rare diseases
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

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== IDENTIFIERS ========== -->

      <fn:string key="id">
        <xsl:value-of select="concat('Orphanet:', fn:number[@key='OrphaCode'])"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="fn:number[@key='OrphaCode']">
          <fn:number key="orpha_code">
            <xsl:value-of select="fn:number[@key='OrphaCode']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="orpha_code"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='Name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Definition'])">
          <fn:null key="definition"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="definition">
            <xsl:value-of select="fn:string[@key='Definition']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CLASSIFICATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='OrphaGroup'])">
          <fn:null key="orpha_group"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="orpha_group">
            <xsl:value-of select="fn:string[@key='OrphaGroup']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:map[@key='DisorderType']">
          <xsl:copy-of select="fn:map[@key='DisorderType']"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="disorder_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE ASSOCIATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='AssociationType'])">
          <fn:null key="association_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="association_type">
            <xsl:value-of select="fn:string[@key='AssociationType']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='AssociationStatus'])">
          <fn:null key="association_status"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="association_status">
            <xsl:value-of select="fn:string[@key='AssociationStatus']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SYNONYMS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='Synonym']/fn:string">
          <fn:array key="synonyms">
            <xsl:for-each select="fn:array[@key='Synonym']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="synonyms"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CROSS-REFERENCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='ExternalReference']/fn:string">
          <fn:array key="xref">
            <xsl:for-each select="fn:array[@key='ExternalReference']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="xref"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Orphanet</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='OrphaCode']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
