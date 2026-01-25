<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  DSLD -> Unified Dietary Supplements Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Dietary Supplements Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ dsld_id                 │ product_id      │ Prefix       │
  │ product_name            │ product_name    │ Direct       │
  │ brand_name              │ brand           │ Direct       │
  │ product_type            │ form            │ Direct       │
  │ ingredients             │ ingredients     │ Restructure  │
  │ manufacturer            │ manufacturer    │ Direct       │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - DSLD provides label data, not verified content
  - Ingredient amounts from label claims
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
      <!-- Primary identifier with DSLD prefix -->
      <fn:string key="product_id">
        <xsl:value-of select="concat('DSLD:', fn:number[@key='dsld_id'])"/>
      </fn:string>

      <!-- Product name -->
      <fn:string key="product_name">
        <xsl:value-of select="fn:string[@key='product_name']"/>
      </fn:string>

      <!-- Brand -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='brand_name']">
          <fn:string key="brand">
            <xsl:value-of select="fn:string[@key='brand_name']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="brand"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Product form -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='product_type']">
          <fn:string key="form">
            <xsl:value-of select="fn:string[@key='product_type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="form"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Serving info -->
      <xsl:if test="fn:string[@key='serving_size']">
        <fn:string key="serving_size">
          <xsl:value-of select="fn:string[@key='serving_size']"/>
        </fn:string>
      </xsl:if>

      <!-- UPC -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='upc']">
          <fn:string key="upc">
            <xsl:value-of select="fn:string[@key='upc']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="upc"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Ingredients -->
      <xsl:if test="fn:array[@key='ingredients']">
        <fn:array key="ingredients">
          <xsl:for-each select="fn:array[@key='ingredients']/fn:map">
            <fn:map>
              <fn:string key="name">
                <xsl:value-of select="fn:string[@key='ingredient_name']"/>
              </fn:string>
              <xsl:if test="fn:string[@key='amount']">
                <fn:string key="amount">
                  <xsl:value-of select="fn:string[@key='amount']"/>
                </fn:string>
              </xsl:if>
              <xsl:if test="fn:string[@key='unit']">
                <fn:string key="unit">
                  <xsl:value-of select="fn:string[@key='unit']"/>
                </fn:string>
              </xsl:if>
              <xsl:if test="fn:string[@key='daily_value_percent']">
                <fn:string key="daily_value_percent">
                  <xsl:value-of select="fn:string[@key='daily_value_percent']"/>
                </fn:string>
              </xsl:if>
              <xsl:if test="fn:string[@key='ingredient_group']">
                <fn:string key="category">
                  <xsl:value-of select="fn:string[@key='ingredient_group']"/>
                </fn:string>
              </xsl:if>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Manufacturer -->
      <xsl:if test="fn:map[@key='manufacturer']">
        <fn:map key="manufacturer">
          <fn:string key="name">
            <xsl:value-of select="fn:map[@key='manufacturer']/fn:string[@key='name']"/>
          </fn:string>
          <xsl:if test="fn:map[@key='manufacturer']/fn:string[@key='country']">
            <fn:string key="country">
              <xsl:value-of select="fn:map[@key='manufacturer']/fn:string[@key='country']"/>
            </fn:string>
          </xsl:if>
        </fn:map>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">DSLD</fn:string>
        <fn:string key="version">v8</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:number[@key='dsld_id']"/>
        </fn:string>
        <fn:string key="data_type">label</fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>
