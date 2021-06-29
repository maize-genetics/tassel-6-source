import org.gradle.api.JavaVersion.VERSION_11
import org.jetbrains.dokka.gradle.DokkaTask
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

// Note Kotlin version needs to be updated in both the buildscript and plugins.
// Dependencies will follow the buildscript

group = "net.maizegenetics"
version = "6.0"

buildscript {
    val kotlinVersion by extra ("1.4.32")

    repositories {
        mavenCentral()
        gradlePluginPortal()
        maven("https://dl.bintray.com/kotlin/kotlin-eap")
    }

    dependencies {
        classpath("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")
        classpath(kotlin("serialization", version = kotlinVersion))
    }
}


plugins {

    val kotlinVersion = "1.4.32"
    java
    kotlin("jvm") version kotlinVersion
    kotlin("plugin.serialization") version kotlinVersion

    // Shadow allows for the creation of fat jars (all dependencies)
    id("com.github.johnrengelman.shadow") version "5.2.0"
    id("application")
    id("org.openjfx.javafxplugin") version "0.0.9"

    id("org.jetbrains.dokka") version "0.10.1"
    `java-library`
    `maven-publish`
    signing

}

application {
    mainClassName = "net.maizegenetics.pipeline.TasselPipeline"
}

apply {
    plugin("kotlinx-serialization")
}

repositories {
    mavenCentral()
    jcenter()
    maven("https://maven.imagej.net/content/groups/public/")
    maven("https://jitpack.io")
    maven("https://dl.bintray.com/kotlin/kotlin-eap")
    maven("https://kotlin.bintray.com/kotlinx")
    maven("https://oss.sonatype.org/content/repositories/snapshots/")
}

javafx {
    modules("javafx.controls", "javafx.web")
}

dependencies {

    val kotlinVersion = rootProject.extra["kotlinVersion"]

    implementation("log4j:log4j:1.2.17")
    implementation("com.google.guava:guava:29.0-jre")
    implementation("com.github.samtools:htsjdk:2.24.1")
    implementation("com.googlecode.efficient-java-matrix-library:ejml:0.23")
    implementation("colt:colt:1.2.0")
    implementation("com.googlecode.json-simple:json-simple:1.1.1")
    implementation("org.glassfish:javax.json:1.1.4")
    implementation("commons-codec:commons-codec:1.11")
    implementation("org.apache.commons:commons-math3:3.6.1")

    implementation("org.openjfx:javafx-controls:11.0.2")
    implementation("org.openjfx:javafx-web:11.0.2")

    implementation("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.4.3")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-core:1.1.0")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-json:1.1.0")

    implementation("khttp:khttp:1.0.0")

    val kotestVersion = "4.2.6"
    listOf("runner-junit5", "assertions-core", "property").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }

}

java {
    sourceCompatibility = VERSION_11
    targetCompatibility = VERSION_11
    withSourcesJar()
}

tasks.withType<KotlinCompile>().configureEach {
    kotlinOptions.jvmTarget = "11"
}

tasks {
    println("Source directories: ${sourceSets["main"].allSource.srcDirs}")
    val dokka by getting(DokkaTask::class) {
        //outputFormat = "html"
        outputFormat = "gfm"
        outputDirectory = "$buildDir/dokka"
        configuration {
            includes = listOf("/Users/edbuckler/Code/biokotlin/src/main/kotlin/biokotlin/packages.md")
        }
    }
}

val dokkaJavadoc = tasks.register<DokkaTask>("dokkaJavadoc") {
    outputFormat = "javadoc"
    outputDirectory = "$buildDir/dokkaJavadoc"
}

val dokkaJar by tasks.creating(Jar::class) {
    group = JavaBasePlugin.DOCUMENTATION_GROUP
    description = "BioKotlin: ${property("version")}"
    archiveClassifier.set("javadoc")
    from(dokkaJavadoc)
}

tasks.test {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
}

publishing {
    publications {

        create<MavenPublication>("maven") {
            artifactId = "biokotlin"
            description = "org.biokotlin:biokotlin:$version"

            from(components["java"])
            artifact(dokkaJar)

            versionMapping {
                usage("java-api") {
                    fromResolutionOf("runtimeClasspath")
                }
                usage("java-runtime") {
                    fromResolutionResult()
                }
            }

            repositories {
                maven {
                    val releasesRepoUrl = "https://oss.sonatype.org/service/local/staging/deploy/maven2"
                    val snapshotsRepoUrl = "https://oss.sonatype.org/content/repositories/snapshots"
                    url = uri(if (version.toString().endsWith("SNAPSHOT")) snapshotsRepoUrl else releasesRepoUrl)
                    credentials {
                        try {
                            username = property("ossrhUsername") as String?
                            password = property("ossrhPassword") as String?
                        } catch (e: Exception) {
                            println("Unable to get username and password for nexus maven central release.")
                        }
                    }
                }
            }

            pom {
                name.set("BioKotlin")
                description.set("BioKotlin aims to be a high-performance bioinformatics library that brings the power and speed of compiled programming languages to scripting and big data environments.")
                url.set("http://www.biokotlin.org/")
                licenses {
                    license {
                        name.set("The Apache License, Version 2.0")
                        url.set("http://www.apache.org/licenses/LICENSE-2.0.txt")
                    }
                }
                developers {
                    developer {
                        name.set("Ed Buckler")
                        email.set("esb33@cornell.edu")
                    }
                    developer {
                        name.set("Vaishnavi Gupta")
                        email.set("vg222@cornell.edu")
                    }
                    developer {
                        id.set("tmc46")
                        name.set("Terry Casstevens")
                        email.set("tmc46@cornell.edu")
                    }
                    developer {
                        name.set("Zack Miller")
                        email.set("zrm22@cornell.edu")
                    }
                    developer {
                        name.set("Lynn Johnson")
                        email.set("lcj34@cornell.edu")
                    }
                    developer {
                        name.set("Brandon Monier")
                        email.set("bm646@cornell.edu")
                    }
                    developer {
                        name.set("Peter Bradbury")
                        email.set("pjb39@cornell.edu")
                    }
                }
                scm {
                    connection.set("scm:git:git://bitbucket.org:bucklerlab/biokotlin.git")
                    developerConnection.set("scm:git:ssh://bitbucket.org:bucklerlab/biokotlin.git")
                    url.set("https://bitbucket.org/bucklerlab/biokotlin/src")
                }
            }
        }
    }
}

signing {
    useGpgCmd()
    sign(publishing.publications["maven"])
}

tasks.javadoc {
    dependsOn("dokkaJavadoc")
    if (JavaVersion.current().isJava9Compatible) {
        (options as StandardJavadocDocletOptions).addBooleanOption("html5", true)
    }
}

