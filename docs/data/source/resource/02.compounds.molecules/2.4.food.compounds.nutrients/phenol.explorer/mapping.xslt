<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Phenol-Explorer -> Food Compounds Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Phenol-Explorer Polyphenol Database)
  Target: ../schema.json (2.4 Food Compounds and Nutrients Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ compound_id                 │ compound_id         │ Direct       │
  │ name                        │ name                │ Direct       │
  │ smiles                      │ smiles              │ Direct       │
  │ molecular_weight            │ molecular_weight    │ Direct       │
  │ class                       │ compound_class      │ Direct       │
  │ subclass                    │ subclass            │ Direct       │
  │ pubchem_cid                 │ pubchem_cid         │ Direct       │
  │ cas_number                  │ cas_number          │ Direct       │
  │ metabolites                 │ metabolites         │ Object array │
  │ pharmacokinetics            │ pharmacokinetics    │ Object       │
  │ food_content                │ food_sources        │ Object array │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - Phenol-Explorer focuses on dietary polyphenols
  - Includes metabolite data (Phase I/II, microbial)
  - Pharmacokinetics includes Cmax, Tmax, AUC
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

      <fn:string key="compound_id">
        <xsl:value-of select="fn:string[@key='compound_id']"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="fn:number[@key='pubchem_cid']">
          <fn:number key="pubchem_cid">
            <xsl:value-of select="fn:number[@key='pubchem_cid']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pubchem_cid"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cas_number'])">
          <fn:null key="cas_number"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cas_number">
            <xsl:value-of select="fn:string[@key='cas_number']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="not(fn:number[@key='molecular_weight']) or fn:number[@key='molecular_weight'] &lt;= 0">
          <fn:null key="molecular_weight"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='molecular_weight']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== POLYPHENOL CLASSIFICATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='class'])">
          <fn:null key="compound_class"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="compound_class">
            <xsl:value-of select="fn:string[@key='class']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='subclass'])">
          <fn:null key="subclass"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="subclass">
            <xsl:value-of select="fn:string[@key='subclass']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METABOLITES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='metabolites']/fn:map">
          <fn:array key="metabolites">
            <xsl:for-each select="fn:array[@key='metabolites']/fn:map">
              <fn:map>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='name']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='type']">
                  <fn:string key="type">
                    <xsl:value-of select="fn:string[@key='type']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='phase']">
                  <fn:string key="phase">
                    <xsl:value-of select="fn:string[@key='phase']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="metabolites"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PHARMACOKINETICS ========== -->

      <xsl:choose>
        <xsl:when test="fn:map[@key='pharmacokinetics']">
          <fn:map key="pharmacokinetics">
            <xsl:if test="fn:map[@key='pharmacokinetics']/fn:number[@key='cmax']">
              <fn:number key="cmax">
                <xsl:value-of select="fn:map[@key='pharmacokinetics']/fn:number[@key='cmax']"/>
              </fn:number>
            </xsl:if>
            <xsl:if test="fn:map[@key='pharmacokinetics']/fn:number[@key='tmax']">
              <fn:number key="tmax">
                <xsl:value-of select="fn:map[@key='pharmacokinetics']/fn:number[@key='tmax']"/>
              </fn:number>
            </xsl:if>
            <xsl:if test="fn:map[@key='pharmacokinetics']/fn:number[@key='auc']">
              <fn:number key="auc">
                <xsl:value-of select="fn:map[@key='pharmacokinetics']/fn:number[@key='auc']"/>
              </fn:number>
            </xsl:if>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pharmacokinetics"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== FOOD SOURCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='food_content']/fn:map">
          <fn:array key="food_sources">
            <xsl:for-each select="fn:array[@key='food_content']/fn:map">
              <fn:map>
                <fn:string key="food_name">
                  <xsl:value-of select="fn:string[@key='food_name']"/>
                </fn:string>
                <xsl:if test="fn:number[@key='content_value']">
                  <fn:number key="content_value">
                    <xsl:value-of select="fn:number[@key='content_value']"/>
                  </fn:number>
                </xsl:if>
                <xsl:if test="fn:string[@key='content_unit']">
                  <fn:string key="content_unit">
                    <xsl:value-of select="fn:string[@key='content_unit']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="food_sources"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Phenol-Explorer</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='compound_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
