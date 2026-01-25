<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  HMDB -> Unified Metabolomics Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Metabolomics Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ accession               │ metabolite_id   │ Direct       │
  │ name                    │ metabolite_name │ Direct       │
  │ chemical_formula        │ formula         │ Direct       │
  │ smiles                  │ smiles          │ Direct       │
  │ inchikey                │ inchikey        │ Direct       │
  │ concentrations          │ concentrations  │ Restructure  │
  │ diseases                │ disease_assoc   │ Restructure  │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - HMDB is the most comprehensive human metabolome database
  - Maps directly to unified metabolomics schema
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
      <!-- Primary identifier -->
      <fn:string key="metabolite_id">
        <xsl:value-of select="fn:string[@key='accession']"/>
      </fn:string>

      <!-- Metabolite name -->
      <fn:string key="metabolite_name">
        <xsl:value-of select="fn:string[@key='name']"/>
      </fn:string>

      <!-- Chemical properties -->
      <xsl:if test="fn:string[@key='chemical_formula']">
        <fn:string key="formula">
          <xsl:value-of select="fn:string[@key='chemical_formula']"/>
        </fn:string>
      </xsl:if>

      <xsl:if test="fn:number[@key='average_molecular_weight']">
        <fn:number key="molecular_weight">
          <xsl:value-of select="fn:number[@key='average_molecular_weight']"/>
        </fn:number>
      </xsl:if>

      <xsl:if test="fn:string[@key='smiles']">
        <fn:string key="smiles">
          <xsl:value-of select="fn:string[@key='smiles']"/>
        </fn:string>
      </xsl:if>

      <xsl:if test="fn:string[@key='inchikey']">
        <fn:string key="inchikey">
          <xsl:value-of select="fn:string[@key='inchikey']"/>
        </fn:string>
      </xsl:if>

      <xsl:if test="fn:string[@key='cas_registry_number']">
        <fn:string key="cas_number">
          <xsl:value-of select="fn:string[@key='cas_registry_number']"/>
        </fn:string>
      </xsl:if>

      <!-- Taxonomy -->
      <xsl:if test="fn:map[@key='taxonomy']">
        <fn:map key="classification">
          <xsl:if test="fn:map[@key='taxonomy']/fn:string[@key='superclass']">
            <fn:string key="superclass">
              <xsl:value-of select="fn:map[@key='taxonomy']/fn:string[@key='superclass']"/>
            </fn:string>
          </xsl:if>
          <xsl:if test="fn:map[@key='taxonomy']/fn:string[@key='class']">
            <fn:string key="class">
              <xsl:value-of select="fn:map[@key='taxonomy']/fn:string[@key='class']"/>
            </fn:string>
          </xsl:if>
        </fn:map>
      </xsl:if>

      <!-- Biofluid locations -->
      <xsl:if test="fn:map[@key='biological_properties']/fn:array[@key='biofluid_locations']">
        <fn:array key="biofluids">
          <xsl:for-each select="fn:map[@key='biological_properties']/fn:array[@key='biofluid_locations']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Concentrations -->
      <xsl:if test="fn:array[@key='concentrations']">
        <fn:array key="concentrations">
          <xsl:for-each select="fn:array[@key='concentrations']/fn:map">
            <fn:map>
              <fn:string key="biofluid">
                <xsl:value-of select="fn:string[@key='biofluid']"/>
              </fn:string>
              <xsl:if test="fn:number[@key='concentration_value']">
                <fn:number key="value">
                  <xsl:value-of select="fn:number[@key='concentration_value']"/>
                </fn:number>
              </xsl:if>
              <fn:string key="unit">
                <xsl:value-of select="fn:string[@key='concentration_units']"/>
              </fn:string>
              <fn:string key="condition">
                <xsl:value-of select="fn:string[@key='condition']"/>
              </fn:string>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Disease associations -->
      <xsl:if test="fn:array[@key='diseases']">
        <fn:array key="disease_associations">
          <xsl:for-each select="fn:array[@key='diseases']/fn:map">
            <fn:map>
              <fn:string key="disease">
                <xsl:value-of select="fn:string[@key='disease_name']"/>
              </fn:string>
              <xsl:if test="fn:string[@key='omim_id']">
                <fn:string key="omim_id">
                  <xsl:value-of select="fn:string[@key='omim_id']"/>
                </fn:string>
              </xsl:if>
              <fn:string key="change">
                <xsl:value-of select="fn:string[@key='concentration_change']"/>
              </fn:string>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">HMDB</fn:string>
        <fn:string key="version">5.0</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='accession']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
